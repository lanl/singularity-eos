//------------------------------------------------------------------------------
// © 2021-2023. Triad National Security, LLC. All rights reserved.  This
// program was produced under U.S. Government contract 89233218CNA000001
// for Los Alamos National Laboratory (LANL), which is operated by Triad
// National Security, LLC for the U.S.  Department of Energy/National
// Nuclear Security Administration. All rights in the program are
// reserved by Triad National Security, LLC, and the U.S. Department of
// Energy/National Nuclear Security Administration. The Government is
// granted for itself and others acting on its behalf a nonexclusive,
// paid-up, irrevocable worldwide license in this material to reproduce,
// prepare derivative works, distribute copies to the public, perform
// publicly and display publicly, and to permit others to do so.
//------------------------------------------------------------------------------

/*
  Compare SpinerEOS to EOSPAC
  on output data text files
  Author: Jonah Miller (jonahm@lanl.gov)
 */

#ifdef SPINER_USE_HDF
#ifdef SINGULARITY_USE_EOSPAC

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <hdf5.h>
#include <hdf5_hl.h>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>

#include <spiner/databox.hpp>
#include <spiner/interpolation.hpp>

#include <singularity-eos/eos/eos.hpp>
#include <singularity-eos/eos/eos_builder.hpp>

#include <eospac-wrapper/eospac_wrapper.hpp>

using namespace singularity;
using namespace EospacWrapper;

using duration = std::chrono::duration<long double>;
using dvec = std::vector<double>;
using ivec = std::vector<int>;
using Spiner::DataBox;
using Spiner::getOnDeviceDataBox;
using Spiner::RegularGrid1D;

#ifdef PORTABILITY_STRATEGY_KOKKOS
using RView = Kokkos::View<Real *>;
using RMirror = typename RView::HostMirror;
#endif

struct Data {
  DataBox rhos, sies, vfracs, Ts;
  bool onDevice = false;
  Data GetOnDevice() const {
    Data other;
    other.rhos = getOnDeviceDataBox(rhos);
    other.sies = getOnDeviceDataBox(sies);
    other.vfracs = getOnDeviceDataBox(vfracs);
    other.Ts = getOnDeviceDataBox(Ts);
    other.onDevice = true;
  }
  void Finalize() {
#ifdef PORTABILITY_STRATEGY_KOKKOS
    using HS = Kokkos::HostSpace;
    using DMS = Kokkos::DefaultExecutionSpace::memory_space;
    constexpr const bool execution_is_host{std::is_same<DMS, HS>::value};
    if (!execution_is_host && onDevice) {
      PORTABLE_FREE(rhos.data());
      PORTABLE_FREE(sies.data());
      PORTABLE_FREE(vfracs.data());
      PORTABLE_FREE(Ts.data());
    }
#endif
  }
};

inline auto now() { return std::chrono::high_resolution_clock::now(); }

template <typename Function>
inline auto get_duration(Function function) {
  auto start = now();
  function();
  auto stop = now();
  return std::chrono::duration_cast<duration>(stop - start).count();
}

void load_data(const std::string &filename, Data &d);
void filter_data_for_EOSPAC(const Data &d_in, std::vector<std::vector<EOS_REAL>> &rhos,
                            std::vector<std::vector<EOS_REAL>> &sies,
                            std::vector<std::vector<EOS_REAL>> &Ts);
void load_eospac_tables(const std::vector<int> &matids,
                        std::vector<EOS_INTEGER> &table_handles);
void load_eos_spiner(const std::string &filename, const std::vector<int> &matids,
                     std::vector<EOS> &eos_h, std::vector<EOS> &eos_d);
void profile_eospac(int ncycles, const std::vector<int> &chunk_sizes, const Data &data,
                    std::vector<EOS_INTEGER> &table_handles, std::vector<double> &timings,
                    std::vector<double> &throughputs, std::vector<Real> &errors);
void profile_spiner_host(int ncycles, const std::vector<int> &chunk_sizes, Data &data,
                         std::vector<EOS> &eos, std::vector<double> &timings,
                         std::vector<double> &throughputs, std::vector<Real> &errors);
void compare_eospac_to_spiner(Data &data, int m, int nside, const std::string &outname,
                              std::vector<EOS_INTEGER> &table_handle,
                              std::vector<EOS> &eos);
void make_table_sim_data(std::ostream &out, const std::string &dumpname, int ncycles,
                         const std::vector<int> &chunk_sizes,
                         const std::vector<Real> &errors_eospac,
                         const std::vector<Real> &errors_spiner,
                         const std::vector<double> &timings_eospac,
                         const std::vector<double> &timings_spiner,
                         const std::vector<double> &throughputs_eospac,
                         const std::vector<double> &throughputs_spiner);
void compare_on_sim_data(int ncycles, const std::vector<int> &chunk_sizes);

void load_data(const std::string &filename, Data &d) {

  std::vector<std::vector<Real>> matrix;

  std::ifstream file;
  std::string line;

  file.open(filename, std::ios::in);
  {
    std::getline(file, line); // first line
    while (std::getline(file, line)) {
      std::stringstream ss(line);
      std::vector<Real> row;
      Real val;
      while (ss >> val) {
        row.push_back(val);
      }
      matrix.push_back(row);
    }
  }
  file.close();

  d.Ts.resize(matrix.size());
  for (int i = 0; i < matrix.size(); i++) {
    int j = matrix[i].size() - 4;
    d.Ts(i) = matrix[i][j];
  }

  int stride = 3;
  int nmats = (matrix[0].size() - 4) / stride;
  d.rhos.resize(nmats, matrix.size());
  d.sies.resize(nmats, matrix.size());
  d.vfracs.resize(nmats, matrix.size());
  for (int i = 0; i < matrix.size(); i++) {
    for (int mat = 0; mat < nmats; mat++) {
      Real frac_mass = matrix[i][mat * stride];
      Real frac_vol = matrix[i][mat * stride + 1];
      Real frac_eng = matrix[i][mat * stride + 2];
      d.rhos(mat, i) = frac_mass / frac_vol;
      d.sies(mat, i) = frac_eng / frac_mass;
      d.vfracs(mat, i) = frac_vol;
    }
  }

  return;
}

// this one must be ragged
void filter_data_for_EOSPAC(const Data &d_in, std::vector<std::vector<EOS_REAL>> &rhos,
                            std::vector<std::vector<EOS_REAL>> &sies,
                            std::vector<std::vector<EOS_REAL>> &Ts) {
  int nmats = d_in.rhos.dim(2);
  int ncells_tot = d_in.rhos.dim(1);
  rhos.resize(nmats);
  sies.resize(nmats);
  Ts.resize(nmats);
  for (int mat = 0; mat < nmats; mat++) {
    int cell_count = 0;
    for (int i = 0; i < ncells_tot; i++) {
      cell_count += (d_in.vfracs(mat, i) > 0) ? 1 : 0;
    }
    rhos[mat].resize(cell_count);
    sies[mat].resize(cell_count);
    Ts[mat].resize(cell_count);
    int c = 0;
    for (int i = 0; i < ncells_tot; i++) {
      if (d_in.vfracs(mat, i) > 0) {
        rhos[mat][c] = densityToSesame(d_in.rhos(mat, i));
        sies[mat][c] = sieToSesame(d_in.sies(mat, i));
        Ts[mat][c] = temperatureToSesame(d_in.Ts(i));
        c++;
      }
    }
  }
}

void load_eospac_tables(const std::vector<int> &matids,
                        std::vector<EOS_INTEGER> &table_handles) {
  constexpr EOS_INTEGER NTABLES = 1;
  EOS_INTEGER tableType[] = {EOS_T_DUt};

  table_handles.resize(matids.size());
  for (int i = 0; i < matids.size(); i++) {
    eosSafeLoad(NTABLES, matids[i], tableType, &(table_handles[i]), true);
  }
}

void load_eos_spiner(const std::string &filename, const std::vector<int> &matids,
                     std::vector<EOS> &eos_h, std::vector<EOS> &eos_d) {
  eos_h.clear();
  eos_d.clear();
  EOSBuilder::EOSType type = EOSBuilder::EOSType::SpinerEOSDependsRhoT;
  EOSBuilder::params_t params;
  params["filename"] = filename;
  params["reproducibility_mode"] = false;

  for (int matid : matids) {
    if (matid > 0) { // non-analytic
      params["matid"] = matid;
      eos_h.push_back(EOSBuilder::buildEOS(type, params));
      eos_d.push_back(eos_h.back().GetOnDevice());
    }
  }
}

void profile_eospac(int ncycles, const std::vector<int> &chunk_sizes, const Data &data,
                    std::vector<EOS_INTEGER> &table_handles, std::vector<double> &timings,
                    std::vector<double> &throughputs, std::vector<Real> &errors) {

  timings.clear();
  errors.clear();
  throughputs.clear();

  int nmats = table_handles.size();
  int npoints_tot = 0;
  std::vector<std::vector<Real>> Ts_out(nmats);
  std::vector<std::vector<EOS_REAL>> Ts_true;

  for (int chunk_size : chunk_sizes) {
    std::cout << "\t\t...chunk size = " << chunk_size << std::endl;
    auto time = get_duration([&]() {
      for (int cycle = 0; cycle < ncycles; cycle++) {

        std::vector<std::vector<EOS_REAL>> rhos;
        std::vector<std::vector<EOS_REAL>> sies;
        std::vector<std::vector<Real>> dx(nmats), dy(nmats);
        filter_data_for_EOSPAC(data, rhos, sies, Ts_true);
        npoints_tot = 0;
        for (int m = 0; m < nmats; m++) {
          int npoints_mat = rhos[m].size();
          dx[m].resize(npoints_mat);
          dy[m].resize(npoints_mat);
          Ts_out[m].resize(npoints_mat);
          npoints_tot += npoints_mat;
        }

        for (int m = 0; m < nmats; m++) {
          int size = rhos[m].size();
          int stride = std::min(size, chunk_size);
          int ncalls = rhos[m].size() / stride;
          EOS_INTEGER nXYPairs = stride;
          for (int c = 0; c < ncalls; c++) {
            eosSafeInterpolate(&(table_handles[m]), nXYPairs, &(rhos[m][c * stride]),
                               &(sies[m][c * stride]), &(Ts_out[m][c * stride]),
                               &(dx[m][c * stride]), &(dy[m][c * stride]), "TofRE",
                               false);
          }
          int n_leftover = rhos[m].size() - stride * ncalls;
          if (n_leftover > 0) {
            EOS_INTEGER nXYPairs = n_leftover;
            int c = stride * ncalls;
            eosSafeInterpolate(&(table_handles[m]), nXYPairs, &(rhos[m][c]),
                               &(sies[m][c]), &(Ts_out[m][c]), &(dx[m][c]), &(dy[m][c]),
                               "TofRE", false);
          }
        }
      }
    });
    time /= (npoints_tot * ncycles);
    timings.push_back(time);

    auto throughput = 1. / time;
    throughputs.push_back(throughput);

    Real L2_error = 0;
    int npoints = 0;
    for (int m = 0; m < nmats; m++) {
      for (int i = 0; i < Ts_true[m].size(); i++) {
        Real diff = Ts_out[m][i] - Ts_true[m][i];
        Real denom = 0.5 * (Ts_out[m][i] + Ts_true[m][i]);
        L2_error += diff * diff / (denom * denom + 1e-20);
        npoints++;
      }
    }
    L2_error = std::sqrt(L2_error / npoints);
    errors.push_back(L2_error);
  }
}

void profile_spiner_host(int ncycles, const std::vector<int> &chunk_sizes, Data &data,
                         std::vector<EOS> &eos, std::vector<double> &timings,
                         std::vector<double> &throughputs, std::vector<Real> &errors) {
  timings.clear();
  errors.clear();
  throughputs.clear();

  int nmats = std::min((int)eos.size(), (int)data.rhos.dim(2));
  int ncells = data.rhos.dim(1);

  std::vector<Real> lambda_v(2);
  Real *lambda = lambda_v.data();

  DataBox Ts(nmats, ncells);

  for (int chunk_size : chunk_sizes) {
    std::cout << "\t\t...chunk size = " << chunk_size << std::endl;
    int stride = std::min(ncells, chunk_size);
    int ncalls = ncells / stride;
    auto time = get_duration([&]() {
      for (int cycle = 0; cycle < ncycles; cycle++) {
        for (int m = 0; m < nmats; m++) {
          for (int c = 0; c < ncalls; c++) {
            int istart = c * stride;
#pragma omp simd
            for (int i = 0; i < stride; i++) {
              Real rho = data.rhos(m, istart + i);
              Real sie = data.sies(m, istart + i);
              if (data.vfracs(m, istart + i) > 0) {
                Ts(m, istart + i) =
                    eos[m].TemperatureFromDensityInternalEnergy(rho, sie, lambda);
              }
            }
          }
          int n_leftover = ncells - stride * ncalls;
          if (n_leftover > 0) {
            int istart = stride * ncalls;
#pragma omp simd
            for (int i = 0; i < n_leftover; i++) {
              Real rho = data.rhos(m, istart + i);
              Real sie = data.sies(m, istart + i);
              if (data.vfracs(m, istart + i) > 0) {
                Ts(m, istart + i) =
                    eos[m].TemperatureFromDensityInternalEnergy(rho, sie, lambda);
              }
            }
          }
        }
      }
    });
    time /= (ncells * ncycles);
    timings.push_back(time);

    auto throughput = 1. / time;
    throughputs.push_back(throughput);

    // debug
    Real L2_error = 0;
    int npoints = 0;
    for (int m = 0; m < nmats; m++) {
      for (int i = 0; i < ncells; i++) {
        // std::cout << "m,i = " << m << ", " << i << std::endl; // debug
        // std::cout << "vfrac = " << data.vfracs(m,i) << std::endl;
        if (data.vfracs(m, i) > 0) {
          Real diff = Ts(m, i) - data.Ts(i);
          Real denom = 0.5 * (Ts(m, i) + data.Ts(i));
          L2_error += diff * diff / (denom * denom + 1e-20);
          npoints++;
        }
      }
    }
    L2_error = std::sqrt(L2_error / npoints);
    errors.push_back(L2_error);
  }
}

void compare_eospac_to_spiner(Data &data, int m, int nside, const std::string &outname,
                              std::vector<EOS_INTEGER> &table_handle,
                              std::vector<EOS> &eos) {

  std::cout << "\t\t...preparing data for EOSPAC" << std::endl;
  std::vector<std::vector<EOS_REAL>> rhos, sies, scratch;

  filter_data_for_EOSPAC(data, rhos, sies, scratch);

  int ncells = rhos[m].size();
  Real rhoMin = densityFromSesame(*std::min_element(rhos[m].begin(), rhos[m].end()));
  Real rhoMax = densityFromSesame(*std::max_element(rhos[m].begin(), rhos[m].end()));
  Real sieMin = sieFromSesame(*std::min_element(sies[m].begin(), sies[m].end()));
  Real sieMax = sieFromSesame(*std::max_element(sies[m].begin(), sies[m].end()));
  std::cout << "\t\t\t...ncells = " << ncells << std::endl;
  std::cout << "\t\t\t...rho = [ " << rhoMin << ", " << rhoMax << "]" << std::endl;
  std::cout << "\t\t\t...sie = [ " << sieMin << ", " << sieMax << "]" << std::endl;

  std::vector<Real> dx(ncells);
  std::vector<Real> dy(ncells);
  std::vector<Real> Ts(ncells);

  std::cout << "\t\t...preparing histogram" << std::endl;
  DataBox EOSPAC((nside + 1), (nside - 1));
  DataBox Spiner((nside + 1), (nside - 1));
  DataBox diffs((nside + 1), (nside - 1));
  DataBox counts((nside + 1), (nside - 1));
  EOSPAC.setRange(0, sieMin, sieMax, nside - 1);
  EOSPAC.setRange(1, rhoMin, rhoMax, nside + 1);
  Spiner.copyMetadata(EOSPAC);
  diffs.copyMetadata(EOSPAC);
  counts.copyMetadata(EOSPAC);
  for (int j = 0; j < EOSPAC.dim(2); j++) {
    for (int i = 0; i < EOSPAC.dim(1); i++) {
      counts(j, i) = Spiner(j, i) = EOSPAC(j, i) = 0;
    }
  }

  std::cout << "\t\t...interpolating EOSPAC" << std::endl;
  EOS_INTEGER nXYPairs = ncells;
  eosSafeInterpolate(&(table_handle[m]), nXYPairs, rhos[m].data(), sies[m].data(),
                     Ts.data(), dx.data(), dy.data(), "TofRE", false);

  std::cout << "\t\t...interpolating spiner" << std::endl;
  std::vector<Real> lambda_v(2);
  Real *lambda = lambda_v.data();
  for (int i = 0; i < ncells; i++) {
    double rho = densityFromSesame(rhos[m][i]);
    double sie = sieFromSesame(sies[m][i]);
    int irho = EOSPAC.range(1).index(rho);
    int ie = EOSPAC.range(0).index(sie);
    Real TPAC = temperatureFromSesame(Ts[i]);
    Real TSP = eos[m].TemperatureFromDensityInternalEnergy(rho, sie, lambda);
    Real diff = TPAC - TSP;
    EOSPAC(irho, ie) += TPAC;
    Spiner(irho, ie) += TSP;
    diffs(irho, ie) += diff * diff;
    counts(irho, ie) += 1;
  }

  std::cout << "\t\t...finishing histogram" << std::endl;
  for (int j = 0; j < EOSPAC.dim(2); j++) {
    for (int i = 0; i < EOSPAC.dim(1); i++) {
      EOSPAC(j, i) /= (counts(j, i) + 1e-20);
      Spiner(j, i) /= (counts(j, i) + 1e-20);
      if (diffs(j, i) > 0) {
        diffs(j, i) = std::sqrt(diffs(j, i) / counts(j, i) + 1e-20);
      }
    }
  }

  std::cout << "\t\t...saving to file " << outname << std::endl;
  herr_t status = H5_SUCCESS;
  hid_t file = H5Fcreate(outname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  status += EOSPAC.saveHDF(file, "EOSPAC");
  status += Spiner.saveHDF(file, "Spiner");
  status += diffs.saveHDF(file, "diffs");
  status += H5Fclose(file);
  std::cout << "\t\t...Saved" << std::endl;
}

void make_table_sim_data(std::ostream &out, const std::string &dumpname, int ncycles,
                         const std::vector<int> &chunk_sizes,
                         const std::vector<Real> &errors_eospac,
                         const std::vector<Real> &errors_spiner,
                         const std::vector<double> &timings_eospac,
                         const std::vector<double> &timings_spiner,
                         const std::vector<double> &throughputs_eospac,
                         const std::vector<double> &throughputs_spiner) {
  out << "# input file = " << dumpname << "\n"
      << "# ncyles = " << ncycles << "\n"
      << "# Setup: 1 dump, 1 CPU. Compare EOSPAC to Spiner performance\n"
      << "# [0]: chunk sizes\n"
      << "# [1]: L2-norm error for EOSPAC\n"
      << "# [2]: L2-norm error for Spiner\n"
      << "# [3]: time/point EOSPAC (microseconds/zone-cycle)\n"
      << "# [4]: time/point Spiner (microseconds/zone-cycle)\n"
      << "# [5]: throughput for EOSPAC (zone-cycles/cpu-second)\n"
      << "# [6]: throughput for Spiner (zone-cycles/cpu-second)" << std::endl;
  for (int i = 0; i < chunk_sizes.size(); i++) {
    out << std::setw(20) << chunk_sizes[i] << "\t" << std::setw(20) << errors_eospac[i]
        << "\t" << std::setw(20) << errors_spiner[i] << "\t" << std::setw(20)
        << 1e6 * timings_eospac[i] << "\t" << std::setw(20) << 1e6 * timings_spiner[i]
        << "\t" << std::setw(20) << throughputs_eospac[i] << "\t" << std::setw(20)
        << throughputs_spiner[i] << std::endl;
  }
}

void compare_on_sim_data(int ncycles, const std::vector<int> &chunk_sizes) {
  std::vector<int> matids = {5030, 3337, 4272, 90003, 95501};
  std::string dumpname = "/scratch/users/jonahm/eap/sim_data_dumps/SimData-dmp000605.txt";
  std::string spiner_fname = "../utils/sesame2spiner/sim-data-materials.sp5";
  Data data_h;
  std::vector<EOS_INTEGER> table_handles;
  std::vector<EOS> eos_h, eos_d;

  std::cout << "Comparing on sim data.." << std::endl;

  std::cout << "\t...loading data" << std::endl;
  load_data(dumpname, data_h);

  std::cout << "\t...loading tables" << std::endl;
  load_eospac_tables(matids, table_handles);
  load_eos_spiner(spiner_fname, matids, eos_h, eos_d);
  Data data_d = data_h.GetOnDevice();

  std::cout << "\t...We have " << data_h.rhos.dim(2) << " materials\n"
            << "\t\t...but we only operate on " << eos_h.size() << " of them.\n"
            << "\t...We have " << data_h.rhos.dim(1) << " cells." << std::endl;

  std::cout << "\t...root finding with EOSPAC" << std::endl;
  std::vector<double> timings_eospac, throughputs_eospac;
  std::vector<Real> errors_eospac;
  profile_eospac(ncycles, chunk_sizes, data_h, table_handles, timings_eospac,
                 throughputs_eospac, errors_eospac);

  /*
  std::cout << "#chunk\terrors\ttimings\tthroughputs" << std::endl;
  for (int i = 0; i < chunk_sizes.size(); i++) {
    std::cout << std::setw(20)
              << chunk_sizes[i] << "\t"
              << errors_eospac[i] << "\t"
              << timings_eospac[i] << "\t"
              << throughputs_eospac[i]
              << std::endl;
  }
  */

  std::cout << "\t...root finding with spiner on 1 CPU" << std::endl;
  std::vector<double> timings_host, throughputs_host;
  std::vector<Real> errors_host;
  profile_spiner_host(ncycles, chunk_sizes, data_h, eos_h, timings_host, throughputs_host,
                      errors_host);

  /*
  std::cout << "#chunk\terrors\ttimings\tthroughputs" << std::endl;
  for (int i = 0; i < chunk_sizes.size(); i++) {
    std::cout << std::setw(20)
              << chunk_sizes[i] << "\t"
              << errors_host[i] << "\t"
              << timings_host[i] << "\t"
              << throughputs_host[i]
              << std::endl;
  }
  */
  make_table_sim_data(std::cout, dumpname, ncycles, chunk_sizes, errors_eospac,
                      errors_host, timings_eospac, timings_host, throughputs_eospac,
                      throughputs_host);
  std::ofstream outfile;
  outfile.open("sim_data_timings.dat", std::ios::out | std::ios::trunc);
  make_table_sim_data(outfile, dumpname, ncycles, chunk_sizes, errors_eospac, errors_host,
                      timings_eospac, timings_host, throughputs_eospac, throughputs_host);
  outfile.close();

  data_d.Finalize();
  for (int m = 0; m < eos_d.size(); m++) {
    eos_d[m].Finalize();
  }

  std::cout << "\t...comparing EOSPAC output to Spiner in histogram..." << std::endl;
  compare_eospac_to_spiner(data_h, 2, 256, "sim_data_diffs.sp5", table_handles, eos_h);
}

int main() {

  // compare_on_sim_data(100, {1, 8, 64, 512, 4096, 32768, 262144, 75449});
  compare_on_sim_data(20, {1, 8, 64, 512, 4096, 32768, 262144, 755449});
  // compare_on_sim_data(1, {4096});

  return 0;
}

#endif // EOSPAC
#endif // SPINER_USE_HDF

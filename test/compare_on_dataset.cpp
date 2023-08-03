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
  on Rage output data text files
  Author: Jonah Miller (jonahm@lanl.gov)
  Copyright: LANL
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

#include <eospac-wrapper/eospac_wrapper.hpp>

using namespace singularity;
using namespace EospacWrapper;

using duration = std::chrono::duration<long double>;
using dvec = std::vector<double>;
using ivec = std::vector<int>;
using Spiner::getOnDeviceDataBox;
using Spiner::RegularGrid1D;

using BEOS = SpinerEOSDependsRhoT;

#ifdef PORTABILITY_STRATEGY_KOKKOS
using RView = Kokkos::View<Real *>;
using RMirror = typename RView::HostMirror;
using KokkosAtomic = Kokkos::MemoryTraits<Kokkos::Atomic>;
#endif

const std::vector<int> matids = {
    5762, // He4
    3720, // Aluminum
    7592, // CH cushion
    2024, // Be Tamper
    3070, // Cr inner shell
    1018, // DT
};
const int nmats = matids.size();
const std::string spiner_fname = "../utils/sesame2spiner/dataset-materials.sp5";
const std::string hist_name = "dataset-histograms.sp5";
constexpr int HIST_NSIDE = 500;

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
    return other;
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
                     std::vector<BEOS> &eos_h, std::vector<BEOS> &eos_d);
void getDataSetBounds(const PortableMDArray<Data> &data, RegularGrid1D &rhoBounds,
                      RegularGrid1D &sieBounds);
double profile_eospac(const PortableMDArray<Data> &data,
                      std::vector<EOS_INTEGER> &table_handles, DataBox &hist_EOSPAC);
double profile_spiner_host(const PortableMDArray<Data> &data, PortableMDArray<BEOS> &eos,
                           DataBox &hist);
double profile_spiner_device(const PortableMDArray<Data> &data,
                             PortableMDArray<BEOS> &eos, const int npoints,
                             DataBox &hist);
void save_histograms(DataBox &EOSPAC, DataBox &spiner_h, DataBox &spiner_hm);

inline void set_zero(DataBox &d) {
  for (int i = 0; i < d.size(); i++)
    d(i) = 0;
}
inline void normalize_hist(DataBox &hist, DataBox &counts) {
  for (int m = 0; m < hist.dim(3); m++) {
    for (int j = 0; j < hist.dim(2); j++) {
      for (int i = 0; i < hist.dim(1); i++) {
        if (counts(m, j, i) > 0) {
          hist(m, j, i) /= (counts(m, j, i) + 1e-5);
        }
      }
    }
  }
}

// ======================================================================

int main(int argc, char *argv[]) {

#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::initialize(argc, argv);
#endif
  {
    if (argc <= 1) {
      std::cerr << argv[0] << ": file1 file2 ..." << std::endl;
      return 1;
    }

    std::cout << "Comparing on a set of dumps." << std::endl;

    int nfiles = argc - 1;
    std::vector<Data> data_h_vec(nfiles);
    PortableMDArray<Data> data_h(data_h_vec.data(), nfiles);
#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::View<Data *> data_dv("data", nfiles);
    auto data_hm = Kokkos::create_mirror_view(data_dv);
    PortableMDArray<Data> data_d(data_dv.data(), nfiles);
#else
    PortableMDArray<Data> data_d, data_hm;
    data_d = data_h; // shallow copies of data_h
    data_hm = data_h;
#endif

    std::cout << "Loading files..." << std::endl;
    // must convert to std::string to sort since strings
    // are comparable and char[] are not.
    std::vector<std::string> infiles(nfiles);
    for (int i = 1; i < argc; i++) {
      infiles[i - 1] = argv[i];
    }
    std::sort(infiles.begin(), infiles.end());
    for (int i = 0; i < nfiles; i++) {
      std::cout << "\t..." << argv[i] << std::endl;
      load_data(infiles[i], data_h(i));
    }

    std::cout << "The data contains " << data_h(0).rhos.dim(2) << " materials.\n"
              << "\tWe are using " << matids.size() << " materials." << std::endl;

    std::cout << "Moving data to device..." << std::endl;
    for (int i = 0; i < argc - 1; i++) {
      data_hm(i) = data_h(i).GetOnDevice();
    }
#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::deep_copy(data_dv, data_hm);
#endif

    std::cout << "Loading tables.." << std::endl;
    std::vector<EOS_INTEGER> table_handles;
    std::vector<BEOS> eos_h_vec, eos_hm_vec;
    load_eospac_tables(matids, table_handles);
    load_eos_spiner(spiner_fname, matids, eos_h_vec, eos_hm_vec);
    PortableMDArray<BEOS> eos_h(eos_h_vec.data(), eos_h_vec.size());
#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::View<BEOS *> eos_dv("eos", eos_h_vec.size());
    auto eos_hm = Kokkos::create_mirror_view(eos_dv);
    for (int m = 0; m < eos_hm_vec.size(); m++) {
      eos_hm(m) = eos_hm_vec[m];
    }
    Kokkos::deep_copy(eos_dv, eos_hm);
    PortableMDArray<BEOS> eos_d(eos_dv.data(), eos_h_vec.size());
#else
    PortableMDArray<BEOS> eos_d(eos_h.data(), eos_h_vec.size());
    PortableMDArray<BEOS> eos_hm(eos_hm_vec.data(), eos_hm_vec.size());
#endif

    RegularGrid1D rhoBounds, sieBounds;
    getDataSetBounds(data_h, rhoBounds, sieBounds);
    std::cout << "Bounds are:\n"
              << "\trho: [" << rhoBounds.min() << ", " << rhoBounds.max() << "]\n"
              << "\tsie: [" << sieBounds.min() << ", " << sieBounds.max() << "]"
              << std::endl;

    std::cout << "Profiling EOSPAC..." << std::endl;
    DataBox hist_EOSPAC(nmats, HIST_NSIDE + 1, HIST_NSIDE - 1);
    hist_EOSPAC.setRange(0, sieBounds);
    hist_EOSPAC.setRange(1, rhoBounds);
    double t_EOSPAC = profile_eospac(data_h, table_handles, hist_EOSPAC);

    std::cout << "Profiling spiner on host..." << std::endl;
    DataBox hist_spiner_h(nmats, HIST_NSIDE + 1, HIST_NSIDE - 1);
    hist_spiner_h.setRange(0, sieBounds);
    hist_spiner_h.setRange(1, rhoBounds);
    double t_spiner_h = profile_spiner_host(data_h, eos_h, hist_spiner_h);

    // TODO: implement me
    std::cout << "Profiling spiner on device..." << std::endl;
#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::View<Real *, KokkosAtomic> hist_spiner_dv("hist", nmats * (HIST_NSIDE + 1) *
                                                                  (HIST_NSIDE - 1));
    auto hist_spiner_hmv = Kokkos::create_mirror_view(hist_spiner_dv);
    DataBox hist_spiner_d(hist_spiner_dv.data(), nmats, HIST_NSIDE + 1, HIST_NSIDE - 1);
    DataBox hist_spiner_hm(hist_spiner_hmv.data(), nmats, HIST_NSIDE + 1, HIST_NSIDE - 1);
#else
    DataBox hist_spiner_d(nmats, HIST_NSIDE + 1, HIST_NSIDE - 1);
    DataBox hist_spiner_hm(hist_spiner_d.data(), nmats, HIST_NSIDE + 1, HIST_NSIDE - 1);
#endif
    hist_spiner_d.setRange(0, sieBounds);
    hist_spiner_d.setRange(1, rhoBounds);
    hist_spiner_hm.setRange(0, sieBounds);
    hist_spiner_hm.setRange(1, rhoBounds);
    double t_spiner_d =
        profile_spiner_device(data_d, eos_d, data_h(0).rhos.dim(1), hist_spiner_d);

#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::deep_copy(hist_spiner_hmv, hist_spiner_dv);
#endif

    std::cout << "Saving to file..." << std::endl;
    save_histograms(hist_EOSPAC, hist_spiner_h, hist_spiner_hm);
    std::cout << "Saved." << std::endl;

    std::cout << "\nTimings (microseconds/zone-cycle):\n"
              << "\tEOSPAC        = " << 1e6 * t_EOSPAC << "\n"
              << "\tspiner host   = " << 1e6 * t_spiner_h << "\n"
              << "\tspiner device = " << 1e6 * t_spiner_d << "\n"
              << "\nThroughput (zone-cycles/second):\n"
              << "\tEOSPAC        = " << 1. / t_EOSPAC << "\n"
              << "\tspiner host   = " << 1. / t_spiner_h << "\n"
              << "\tspiner device = " << 1. / t_spiner_d << "\n"
              << std::endl;

    for (int m = 0; m < nmats; m++) {
      eos_hm(m).Finalize();
    }
    for (int i = 0; i < nfiles; i++) {
      data_hm(i).Finalize();
    }
  }
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::finalize();
#endif

  return 0;
}

// ======================================================================

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
      d.rhos(mat, i) = frac_vol > 0 ? frac_mass / (frac_vol + 1e-5) : 0;
      // if (d.rhos(mat,i) < 1e-5) d.rhos(mat,i) = 1e-5;
      d.sies(mat, i) = frac_mass > 0 ? frac_eng / (frac_mass + 1e-5) : 0;
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
                     std::vector<BEOS> &eos_h, std::vector<BEOS> &eos_d) {
  eos_h.clear();
  eos_d.clear();
  for (int matid : matids) {
    if (matid <= 0) continue;
    eos_h.push_back(SpinerEOSDependsRhoT(filename, matid, false));
    eos_d.push_back(eos_h.back().GetOnDevice());
  }
}

void getDataSetBounds(const PortableMDArray<Data> &data, RegularGrid1D &rhoBounds,
                      RegularGrid1D &sieBounds) {
  const int nfiles = data.GetDim1();
  Real rhoMinAll = std::numeric_limits<Real>::infinity();
  Real rhoMaxAll = -std::numeric_limits<Real>::infinity();
  Real sieMinAll = std::numeric_limits<Real>::infinity();
  Real sieMaxAll = -std::numeric_limits<Real>::infinity();
  for (int n = 0; n < nfiles; n++) {
    auto rhoFirst = data(n).rhos.data();
    auto rhoLast = rhoFirst + data(n).rhos.size();
    Real rhoMin = *std::min_element(rhoFirst, rhoLast);
    Real rhoMax = *std::max_element(rhoFirst, rhoLast);
    rhoMinAll = std::min(rhoMin, rhoMinAll);
    rhoMaxAll = std::max(rhoMax, rhoMaxAll);

    auto sieFirst = data(n).sies.data();
    auto sieLast = sieFirst + data(n).sies.size();
    Real sieMin = *std::min_element(sieFirst, sieLast);
    Real sieMax = *std::max_element(sieFirst, sieLast);
    sieMinAll = std::min(sieMin, sieMinAll);
    sieMaxAll = std::max(sieMax, sieMaxAll);
  }
  rhoBounds = RegularGrid1D(rhoMinAll, rhoMaxAll, HIST_NSIDE + 1);
  sieBounds = RegularGrid1D(sieMinAll, sieMaxAll, HIST_NSIDE - 1);
}

double profile_eospac(const PortableMDArray<Data> &data,
                      std::vector<EOS_INTEGER> &table_handles, DataBox &hist_EOSPAC) {

  const int ndumps = data.GetDim1();
  const int nmats = table_handles.size();
  const int npoints = data(0).rhos.dim(1);
  std::vector<std::vector<Real>> Ts_out(nmats), dx(nmats), dy(nmats);
  std::vector<std::vector<EOS_REAL>> rhos, sies, Ts_true;

  DataBox counts;
  counts.copyMetadata(hist_EOSPAC);
  RegularGrid1D rhoBounds = hist_EOSPAC.range(1);
  RegularGrid1D sieBounds = hist_EOSPAC.range(0);

  set_zero(hist_EOSPAC);
  set_zero(counts);

  auto time = get_duration([&]() {
    for (int d = 0; d < ndumps; d++) {
      filter_data_for_EOSPAC(data(d), rhos, sies, Ts_true);
      for (int m = 0; m < nmats; m++) {
        int npoints_mat = rhos[m].size();
        dx[m].resize(npoints_mat);
        dy[m].resize(npoints_mat);
        Ts_out[m].resize(npoints_mat);

        EOS_INTEGER nXYPairs = npoints_mat;
        eosSafeInterpolate(&(table_handles[m]), nXYPairs, &(rhos[m][0]), &(sies[m][0]),
                           &(Ts_out[m][0]), &(dx[m][0]), &(dy[m][0]), "TofRE", false);

        for (int i = 0; i < npoints_mat; i++) {
          Real rho = densityFromSesame(rhos[m][i]);
          Real sie = sieFromSesame(sies[m][i]);
          Real T = temperatureFromSesame(Ts_out[m][i]);
          int iRho = rhoBounds.index(rho);
          int iSie = sieBounds.index(sie);
          hist_EOSPAC(m, iRho, iSie) += T;
          counts(m, iRho, iSie) += 1;
        }
      }
    }
    normalize_hist(hist_EOSPAC, counts);
  });
  return time / (nmats * ndumps * npoints);
}

double profile_spiner_host(const PortableMDArray<Data> &data, PortableMDArray<BEOS> &eos,
                           DataBox &hist) {

  const int ndumps = data.GetDim1();
  const int nmats = eos.GetDim1();
  const int npoints = data(0).rhos.dim(1);

  DataBox counts;
  counts.copyMetadata(hist);
  RegularGrid1D rhoBounds = hist.range(1);
  RegularGrid1D sieBounds = hist.range(0);
  set_zero(hist);
  set_zero(counts);

  DataBox lambda_db(nmats, npoints, 2);

  auto time = get_duration([&]() {
    for (int d = 0; d < ndumps; d++) {
      for (int m = 0; m < nmats; m++) {
#pragma omp simd
        for (int i = 0; i < npoints; i++) {
          Real rho = data(d).rhos(m, i);
          Real sie = data(d).sies(m, i);
          int iRho = rhoBounds.index(rho);
          int iSie = sieBounds.index(sie);
          Real *lambda = &(lambda_db(m, i, 0));
          if (data(d).vfracs(m, i) > 0) {
            Real T = eos(m).TemperatureFromDensityInternalEnergy(rho, sie, lambda);
            hist(m, iRho, iSie) += T;
            counts(m, iRho, iSie) += 1;
          }
        }
      }
    }
    normalize_hist(hist, counts);
  });
  return time / (nmats * ndumps * npoints);
}

double profile_spiner_device(const PortableMDArray<Data> &data,
                             PortableMDArray<BEOS> &eos, const int npoints,
                             DataBox &hist) {
  const int ndumps = data.GetDim1();
  const int nmats = eos.GetDim1();

  RegularGrid1D rhoBounds = hist.range(1);
  RegularGrid1D sieBounds = hist.range(0);

#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::View<Real *, KokkosAtomic> counts_dv("counts", nmats * (HIST_NSIDE + 1) *
                                                             (HIST_NSIDE - 1));
  DataBox counts(counts_dv.data(), nmats, HIST_NSIDE + 1, HIST_NSIDE - 1);
#else
  DataBox counts_h(HIST_NSIDE + 1, HIST_NSIDE - 1);
  DataBox counts = getOnDeviceDataBox(counts_h);
#endif

  DataBox lambda_db_h(nmats, npoints, 2);
  auto lambda_db = getOnDeviceDataBox(lambda_db_h);

  portableFor(
      "zeros", 0, nmats, 0, hist.dim(2), 0, hist.dim(1),
      PORTABLE_LAMBDA(const int k, const int j, const int i) {
        hist(k, j, i) = 0;
        counts(k, j, i) = 0;
      });

  auto time = get_duration([&]() {
    portableFor(
        "spiner", 0, ndumps, 0, nmats, 0, npoints,
        PORTABLE_LAMBDA(const int d, const int m, const int i) {
          Real rho = data(d).rhos(m, i);
          Real sie = data(d).sies(m, i);
          Real vfrac = data(d).vfracs(m, i);
          Real *lambda = &(lambda_db(m, i, 0));
          Real T = eos(m).TemperatureFromDensityInternalEnergy(rho, sie, lambda);
          int iRho = rhoBounds.index(rho);
          int iSie = sieBounds.index(sie);
          if (vfrac > 0) {
            hist(m, iRho, iSie) += T;
            counts(m, iRho, iSie) += 1;
          }
        });
    portableFor(
        "normalize", 0, nmats, 0, hist.dim(2), 0, hist.dim(1),
        PORTABLE_LAMBDA(const int &m, const int &j, const int &i) {
          if (counts(m, j, i) > 0) hist(m, j, i) /= (counts(m, j, i) + 1e-5);
        });
  });
  return time / (nmats * ndumps * npoints);
}

void save_histograms(DataBox &EOSPAC, DataBox &spiner_h, DataBox &spiner_hm) {
  herr_t status = H5_SUCCESS;
  hid_t file = H5Fcreate(hist_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  status += EOSPAC.saveHDF(file, "EOSPAC");
  status += spiner_h.saveHDF(file, "spiner_h");
  status += spiner_hm.saveHDF(file, "spiner_d");
  status += H5Fclose(file);
}

#endif // EOSPAC
#endif // SPINER_USE_HDF

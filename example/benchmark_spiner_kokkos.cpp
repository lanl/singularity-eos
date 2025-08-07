//------------------------------------------------------------------------------
// © 2021-2025. Triad National Security, LLC. All rights reserved.  This
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
  Example usage: ./example/benchmark_spiner_gpu  path/for/results nRho nT ./materials.sp5 [matids]

  Authors: Erin O'Neil and Joshua Basabe
 */

// C++ Headers
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

// Data headers
#include <hdf5.h>
#include <hdf5_hl.h>

// Spiner headers
#include <ports-of-call/portability.hpp>
#include <singularity-eos/base/sp5/singularity_eos_sp5.hpp>
#include <spiner/databox.hpp>
#include <spiner/interpolation.hpp>
#include <spiner/sp5.hpp>
#include <spiner/spiner_types.hpp>

// Get EOS Models
#include "singularity-eos/eos/eos_spiner.hpp"
#include <eospac-wrapper/eospac_wrapper.hpp>
#include <singularity-eos/eos/eos.hpp>

#include "Kokkos_Core.hpp"

using namespace singularity;
using namespace EospacWrapper;
using namespace std::chrono;
#include <experimental/filesystem> //not sure why I had to do this
namespace fs = std::experimental::filesystem;
using DataBox = Spiner::DataBox<Real>;
using RegularGrid1D = Spiner::RegularGrid1D<Real>;
using Real = double;

// == Creates the Bounds of the Grid ==
class Bounds { // check for edge effects
 public:
  Bounds(Real min, Real max, int N) : offset(0.0) {
    constexpr Real epsilon = std::numeric_limits<float>::epsilon();
    const Real min_offset = 10 * std::abs(epsilon);
    if (min <= 1e-16)
      offset = 1.1 * std::abs(min) + min_offset; // changed this from (min <= 0) because
                                                 // small values will make the log blow up
    min += offset;
    max += offset;
    min = std::log10(min);
    max = std::log10(max);
    dx = (max - min) / (N - 1);
    x0 = min;
    size = N;
  }
  Real i2lin(int i) const {
    return std::pow(10.0, x0 + static_cast<Real>(i) * dx) - offset;
  }

 private:
  Real x0, dx, offset;
  int size;
};

// == The following function creates and fills in a csv of grid values ==
void write_matrix_csv(const std::string &filename,
                      const std::vector<std::vector<Real>> &matrix) {
  std::ofstream file(filename);
  file << std::setprecision(14);
  for (const auto &row : matrix) {
    for (size_t j = 0; j < row.size(); ++j) {
      file << row[j];
      if (j < row.size() - 1) file << ",";
    }
    file << "\n";
  }
}

// == The following saves the rho and T values corresponding to the indices (for plotting
// purposes) ==
void write_vector_csv(const std::string &filename, const std::vector<Real> &vec) {
  std::ofstream file(filename);
  file << std::setprecision(14);
  for (const auto &val : vec) {
    file << val << "\n";
  }
}

// == Main loop ===
int main(int argc, char *argv[]) {
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::initialize();
#endif

  if (argc < 5) {
    std::cerr
        << "Usage: " << argv[0]
        << " path/to/output nRho nT sp5_file matid1 matid2 ..."; // Is this good usage?
    return 1;
  }
  fs::path base_output_path(argv[1]);
  // if the folder does not exist, try to make one
  if (!fs::exists(base_output_path)) {
    if (!fs::create_directories(base_output_path)) {
      std::cerr << "Failed to create output directory: " << base_output_path << "\n";
      return 1;
    }
  }
  int nRho = std::stoi(
      argv[2]); // test on batch sizes of 144, 269, 512, 2048 or round powers of ten
  int nT = std::stoi(argv[3]);
  std::string sp5file = argv[4];

  std::vector<int> matids;
  for (int i = 5; i < argc; ++i)
    matids.push_back(std::atoi(argv[i]));

  // == Iterate through each material ==
  for (int matid : matids) {
    const int ntimes = 20; // Perform 20 time trials

    // == Get material metadata ==
    SesameMetadata meta;
    eosGetMetadata(matid, meta);

    // == Set up bounds ==
    Real rhoMin = 1.1 * std::max(meta.rhoMin, 1e-5);
    Real rhoMax = 0.9 * meta.rhoMax;
    Real TMin = 1.1 * std::max(meta.TMin, 1.0);
    Real TMax = 0.9 * meta.TMax;
    Bounds rhoBounds(rhoMin, rhoMax, nRho);
    Bounds TBounds(TMin, TMax, nT);

    SpinerEOSDependsRhoT eos_rt(sp5file, matid);
    SpinerEOSDependsRhoSie eos_rs(
        sp5file,
        matid); // eventually make it SpinerEOSDependsRhoSie<NullTransfom>, for example
    EOSPAC eos_ref(matid);

    // These vectors will store the compute time for each model for the 20 trials
    std::vector<double> time_sie_rt_list;

    Kokkos::View<Real *> rhos_gpu("rhos_d", nRho); // DOES THE FOLLOWING LINES WORK OKAY?
    Kokkos::View<Real *> temps_gpu("temps_d", nT);

    auto rhos_h = Kokkos::create_mirror_view(rhos_gpu);
    auto temps_h = Kokkos::create_mirror_view(temps_gpu);

    for (int i = 0; i < nRho; ++i)
      rhos_h[i] = rhoBounds.i2lin(i);
    for (int j = 0; j < nT; ++j)
      temps_h[j] = TBounds.i2lin(j);

    Kokkos::deep_copy(rhos_gpu, rhos_h);
    Kokkos::deep_copy(temps_gpu, temps_h);

    Kokkos::View<Real **> sie_rt_gpu("sie_rt_d", nT, nRho);

    auto sie_rt_h = Kokkos::create_mirror_view(sie_rt_gpu);

    // move eos to device
    eos_rt = eos_rt.GetOnDevice();

    auto t0 = high_resolution_clock::now(); // e(rho, T)
    Kokkos::parallel_for(
        "sie_rt", Kokkos::MDRangePolicy({0, 0}, {nT, nRho}),
        KOKKOS_LAMBDA(const int it, const int ir) {
          sie_rt_gpu(it, ir) =
              eos_rt.InternalEnergyFromDensityTemperature(rhos_gpu(ir), temps_gpu(it)); //start with only 1 inversion
          //          sie_rt_gpu(ir, it) =
          //              eos_rt.InternalEnergyFromDensityTemperature(rhos_gpu(ir),
          //              temps_gpu(it));
        });
    Kokkos::fence();
    auto t1 = high_resolution_clock::now();
    auto elapsed = t1 - t0;
    time_sie_rt_list.push_back(elapsed.count());
    // may want to time this seperately
    Kokkos::deep_copy(sie_rt_h, sie_rt_gpu);
    // cleanup eos on device
    eos_rt.Finalize();
    // === Save the computed data ===

    // Convert to std::vector
    std::vector<std::vector<Real>> sie_matrix(nT, std::vector<Real>(nRho));
    for (int i = 0; i < nT; ++i)
      for (int j = 0; j < nRho; ++j)
        sie_matrix[i][j] = sie_rt_h(i, j);

    // Build filenames and write files
    std::string sie_filename = (base_output_path / ("sie_rt_matid" + std::to_string(matid) + ".csv")).string();
    write_matrix_csv(sie_filename, sie_matrix);

    // Save rho and T if desired (you could move these outside the loop if the same for all matids)
    std::vector<Real> rho_vals(nRho), T_vals(nT);
    for (int i = 0; i < nRho; ++i) rho_vals[i] = rhos_h(i);
    for (int j = 0; j < nT; ++j) T_vals[j] = temps_h(j);

    write_vector_csv((base_output_path / "rho.csv").string(), rho_vals);
    write_vector_csv((base_output_path / "T.csv").string(), T_vals);

    std::string timing_filename = (base_output_path / ("timing_sie_rt_matid" + std::to_string(matid) + ".csv")).string();
    write_vector_csv(timing_filename, time_sie_rt_list);


    std::cout << "Benchmark complete for material " << std::to_string(matid) << "\n";
  }

  Kokkos::finalize();
  return 0;
}

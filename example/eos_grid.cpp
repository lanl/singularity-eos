//------------------------------------------------------------------------------
// © 2024. Triad National Security, LLC. All rights reserved.  This
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

// C headers
#include <cmath>

// C++ headers
#include <chrono>
#include <iostream>
#include <string>

// HDF5 we'll need this for I/O
#include <hdf5.h>
#include <hdf5_hl.h>

// This library contains portable utilities
#include <ports-of-call/portability.hpp>

// This contains useful tools for preventing things like divide by zero
#include <singularity-eos/base/robust_utils.hpp>
// Needed to import the eos models
#include <singularity-eos/eos/eos.hpp>

// This library contains the spiner table object, which we will use to
// store our output
#include <spiner/databox.hpp>
#include <spiner/interpolation.hpp>

// for timing
using duration = std::chrono::nanoseconds;

// These are the specializations of spiner we will use
using DataBox = Spiner::DataBox<Real>;
using RegularGrid1D = Spiner::RegularGrid1D<Real>;

// Set the EOS you want to use here.
using EOS = singularity::StellarCollapse;

int main(int argc, char *argv[]) {
  // This is needed for Kokkos
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::initialize();
#endif
  // note the scoping here. This means Kokkos objects are cleaned up
  // before finalization
  {
    if (argc < 8) {
      std::cerr << "Usage: " << argv[0]
                << " savefilename log10rhomin log10rhomax numrho log10Tmin log10Tmax "
                   "numT eosargs..."
                << std::endl;
      std::exit(1);
    }
    const std::string savefilename = argv[1];
    const double lrho_min = atof(argv[2]);
    const double lrho_max = atof(argv[3]);
    const int nrho = atoi(argv[4]);
    RegularGrid1D lrhoGrid(lrho_min, lrho_max, nrho);

    const double lT_min = atof(argv[5]);
    const double lT_max = atof(argv[6]);
    const int nT = atoi(argv[7]);
    RegularGrid1D lTGrid(lT_min, lT_max, nT);

    // This is the databox we will evaluate onto
    DataBox press_d(Spiner::AllocationTarget::Device, nrho, nT);
    press_d.setRange(0, lTGrid); // note 0 is rightmost index
    press_d.setRange(1, lrhoGrid);
    // This is the host mirror, we will copy into it
    DataBox press_h(Spiner::AllocationTarget::Host, nrho, nT);
    press_h.setRange(0, lTGrid);
    press_h.setRange(1, lrhoGrid);

    // These are the pre-evaluated density and temperature points
    std::cout << "Filling grid points..." << std::endl;
    DataBox rhos(Spiner::AllocationTarget::Device, nrho);
    portableFor(
        "Set rho", 0, nrho,
        PORTABLE_LAMBDA(const int i) { rhos(i) = std::pow(10., lrhoGrid.x(i)); });
    DataBox Ts(Spiner::AllocationTarget::Device, nT);
    portableFor(
        "Set T", 0, nT,
        PORTABLE_LAMBDA(const int i) { Ts(i) = std::pow(10., lTGrid.x(i)); });

    // if you have arguments you want to pass into your EOS load them here
    const std::string loadname = argv[8];
    const double Ye = atof(argv[9]);

    std::cout << "Initializing EOS..." << std::endl;
    EOS eos_h(loadname);
    EOS eos_d = eos_h.GetOnDevice();

    std::cout << "Evaluating..." << std::endl;
    // start timers
#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::fence();
#endif
    auto start = std::chrono::high_resolution_clock::now();
    portableFor(
        "P(rho, T)", 0, nrho, 0, nT, PORTABLE_LAMBDA(const int j, const int i) {
          // Some EOS objects take lambdas. If you want to set
          // them locally, do it like this. Otherwise, you may
          // want to create a device-side array for them to
          // avoid race conditions.
          Real lambda[2];
          lambda[0] = Ye;
          press_d(j, i) = eos_d.PressureFromDensityTemperature(rhos(j), Ts(i), lambda);
        });
    // stop timers
#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::fence();
#endif
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration_device = std::chrono::duration_cast<duration>(stop - start);

    // copy data from device to host
    std::cout << "Copying to host..." << std::endl;
    portableCopyToHost(press_h.data(), press_d.data(), press_h.sizeBytes());

    // For fun lets also evalaute all on host
    std::cout << "Generate arrays on host..." << std::endl;
    DataBox rhos_h(Spiner::AllocationTarget::Host, nrho);
    for (int i = 0; i < nrho; ++i) {
      rhos_h(i) = std::pow(10., lrhoGrid.x(i));
    }
    DataBox Ts_h(Spiner::AllocationTarget::Host, nT);
    for (int i = 0; i < nT; ++i) {
      Ts_h(i) = std::pow(10., lTGrid.x(i));
    }
    std::cout << "Looping on host" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    for (int j = 0; j < nT; ++j) {
      for (int i = 0; i < nrho; ++i) {
        Real lambda[2];
        lambda[0] = Ye;
        press_h(j, i) = eos_h.PressureFromDensityTemperature(rhos_h(j), Ts_h(i), lambda);
      }
    }
    stop = std::chrono::high_resolution_clock::now();
    auto duration_host = std::chrono::duration_cast<duration>(stop - start);

    std::cout << "Saving file..." << std::endl;
    press_h.saveHDF(savefilename); // if you have HDF5 enabled, thats all it takes
    std::cout << "Time per point, (ns):\n"
              << "host = " << duration_host.count() / (static_cast<Real>(nrho) * nT)
              << "\n"
              << "device = " << duration_device.count() / (static_cast<Real>(nrho) * nT)
              << std::endl;
  }
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::finalize();
#endif
  return 0;
}

//------------------------------------------------------------------------------
// © 2021. Triad National Security, LLC. All rights reserved.  This
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
  Tool for profiling the Stellar Collapse Table Reader
  Authors: Jonah Miller
 */

#ifdef SPINER_USE_HDF

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

// #include "hdf5.h"
// #include "hdf5_hl.h"

#include <fast-math/logs.hpp>
#include <ports-of-call/portability.hpp>
#include <spiner/databox.hpp>
#include <spiner/interpolation.hpp>
#include <spiner/sp5.hpp>
#include <spiner/spiner_types.hpp>

#include <eos/eos.hpp>

using namespace Spiner;
using namespace singularity;

using duration = std::chrono::microseconds;
// constexpr char DIFFS_NAME[] = "diffs.sp5";
constexpr Real MP = 1.67262171e-24; // proton mass
constexpr int NTRIALS = 10;

struct Bounds {
  Bounds() {}
  Bounds(Real min, Real max, int N, Real offset)
      : grid(RegularGrid1D(min, max, N)), offset(offset) {}
  Bounds(Real min, Real max, int N, bool convertToLog = false) : offset(0) {
    if (convertToLog) {
      constexpr Real epsilon = std::numeric_limits<float>::epsilon();
      const Real min_offset =
          std::max(10 * std::abs(epsilon), std::abs(epsilon * max));
      if (min < 0)
        offset = std::abs(min) + min_offset;
      else if (min == 0) {
        offset = min_offset;
      }
      min += offset;
      max += offset;
      min = std::log10(std::abs(min)); // fast, floating-point log
      max = std::log10(std::abs(max));
    }
    grid = RegularGrid1D(min, max, N);
  }
  PORTABLE_INLINE_FUNCTION Real log2lin(Real xl) const {
    return pow(10., xl) - offset;
  }
  PORTABLE_INLINE_FUNCTION Real i2lin(int i) const {
    return log2lin(grid.x(i));
  }
  RegularGrid1D grid;
  Real offset;
};

int main(int argc, char *argv[]) {

  // herr_t status = H5_SUCCESS;

#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::initialize();
#endif
  {
    if (argc < 4) {
      std::cerr << "Usage: " << argv[0] << "filename gamma nfine" << std::endl;
      std::exit(1);
    }
    std::string filename = argv[1];
    Real gamma = std::atof(argv[2]);
    int nfine = std::atoi(argv[3]);

    if (nfine < 1) {
      std::cerr << "We need at least one interpolation point" << std::endl;
      std::exit(1);
    }
    if (gamma <= 1) {
      std::cerr << "gamma - 1 must be a positive number" << std::endl;
      std::exit(1);
    }

    std::cout << "Profiling a stellar collapse table" << std::endl;
    std::cout << "\t...Creating an ideal gas equation of state with gamma = "
              << gamma << " to compare to." << std::endl;

    const Real Cv = 1. / (MP * (gamma - 1));
    IdealGas ig(gamma - 1, Cv);
    auto ig_d = ig.GetOnDevice();

    std::cout << "\t...Loading the stellar collapse table" << std::endl;
    // Don't use SP5 or filter bmod. Filtering bmod improves accuracy but is
    // slow at startup.
    StellarCollapse sc(filename, false, true);
    auto sc_d = sc.GetOnDevice();

    std::cout << "\t...Table bounds are:\n"
              << "\t\tlog(rho) in [" << sc.lRhoMin() << ", " << sc.lRhoMax()
              << "]\n"
              << "\t\tlog(T)   in [" << std::log10(sc.TMin()) << ", "
              << std::log10(sc.TMax()) << "]\n"
              << "\t\tYe       in [" << sc.YeMin() << ", " << sc.YeMax()
              << "]\n"
              << "\t\tsie      in [" << sc.sieMin() << ", " << sc.sieMax()
              << "]" << std::endl;

    std::cout << "\t...Ofsets are:\n"
              << "\t\tlRho = " << sc.lRhoOffset() << "\n"
              << "\t\tlT   = " << sc.lTOffset() << "\n"
              << "\t\tlE   = " << sc.lEOffset() << std::endl;

    std::cout << "\t...Allocating databoxes for diff" << std::endl;
#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::View<Real *> diffs_press_dv("press", nfine * nfine * nfine);
    Kokkos::View<Real *> diffs_t_dv("T", nfine * nfine * nfine);
    auto diffs_press_hv = Kokkos::create_mirror_view(diffs_press_dv);
    auto diffs_t_hv = Kokkos::create_mirror_view(diffs_t_dv);
    DataBox diffs_press_d(diffs_press_dv.data(), nfine, nfine, nfine);
    DataBox diffs_press_h(diffs_press_hv.data(), nfine, nfine, nfine);
    DataBox diffs_t_d(diffs_t_dv.data(), nfine, nfine, nfine);
    DataBox diffs_t_h(diffs_t_hv.data(), nfine, nfine, nfine);
#else  // host only. Add other strategies above as needed
    DataBox diffs_press_d(nfine, nfine, nfine);
    DataBox diffs_t_d(nfine, nfine, nfine);
    DataBox diffs_press_h = diffs_press_d; // shallow copy
    DataBox diffs_t_h = diffs_t_d;
#endif // PORTABILITY_STRATEGY
    DataBox lambdas_d(AllocationTarget::Device, nfine, nfine, nfine, 2);

    std::cout
        << "\t...Setting loop bounds.\n"
        << "\t\tBounds set to force extrapolation and no exact grid points."
        << std::endl;
    Bounds lTBounds(sc.TMin() - 0.1 * std::abs(sc.TMin()),
                    sc.TMax() + 0.1 * sc.TMax(), nfine, true);
    Bounds lRhoBounds(sc.rhoMin() - 0.1 * std::abs(sc.rhoMin()),
                      sc.rhoMax() + 0.1 * sc.rhoMax(), nfine, true);
    Bounds YeBounds(sc.YeMin() - 0.1 * std::abs(sc.YeMin()),
                    sc.YeMax() + 0.1 * sc.YeMax(), nfine, false);
    Bounds lEBounds(sc.sieMin() - 0.1 * std::abs(sc.sieMin()),
                    sc.sieMax() + 0.1 * sc.sieMax(), nfine, true);
    diffs_press_d.setRange(0, lRhoBounds.grid);
    diffs_press_d.setRange(1, lTBounds.grid);
    diffs_press_d.setRange(2, YeBounds.grid);
    diffs_press_h.setRange(0, lRhoBounds.grid);
    diffs_press_h.setRange(1, lTBounds.grid);
    diffs_press_h.setRange(2, YeBounds.grid);

    diffs_t_d.setRange(0, lRhoBounds.grid);
    diffs_t_d.setRange(1, lEBounds.grid);
    diffs_t_d.setRange(2, YeBounds.grid);
    diffs_t_h.setRange(0, lRhoBounds.grid);
    diffs_t_h.setRange(1, lEBounds.grid);
    diffs_t_h.setRange(2, YeBounds.grid);

    std::cout << "\t...Profiling interpolation P(rho, T, Ye)..." << std::endl;
    duration durationInterp = duration::zero();
    for (int trial = 0; trial < NTRIALS; ++trial) { // includes launch latency
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::fence();
#endif
      auto start = std::chrono::high_resolution_clock::now();
      portableFor(
          "pressure from rho, T, Ye", 0, nfine, 0, nfine, 0, nfine,
          PORTABLE_LAMBDA(const int iYe, const int iT, const int irho) {
            const Real rho = lRhoBounds.i2lin(irho);
            const Real T = lTBounds.i2lin(iT);
            const Real Ye = YeBounds.i2lin(iYe);
            lambdas_d(iYe, iT, irho, 0) = Ye;
            diffs_press_d(iYe, iT, irho) = sc_d.PressureFromDensityTemperature(
                rho, T, &(lambdas_d(iYe, iT, irho, 0)));
          });
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::fence();
#endif
      auto stop = std::chrono::high_resolution_clock::now();
      durationInterp += std::chrono::duration_cast<duration>(stop - start);
    }

    std::cout << "\t...Computing difference from ideal gas..." << std::endl;
    portableFor(
        "ideal gas pressure from rho, T, Ye", 0, nfine, 0, nfine, 0, nfine,
        PORTABLE_LAMBDA(const int iYe, const int iT, const int irho) {
          const Real rho = lRhoBounds.i2lin(irho);
          const Real T = lTBounds.i2lin(iT);
          diffs_press_d(iYe, iT, irho) -=
              ig_d.PressureFromDensityTemperature(rho, T);
        });
#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::deep_copy(diffs_press_hv, diffs_press_dv);
#endif

    std::cout << "\t...Profiling root find T(rho, e, Ye)..." << std::endl;
    duration durationRoot = duration::zero();
    for (int trial = 0; trial < NTRIALS; ++trial) { // includes launch latency
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::fence();
#endif
      auto start = std::chrono::high_resolution_clock::now();
      portableFor(
          "Temperature from rho, sie, Ye", 0, nfine, 0, nfine, 0, nfine,
          PORTABLE_LAMBDA(const int iYe, const int iE, const int irho) {
            const Real rho = lRhoBounds.i2lin(irho);
            const Real sie = lEBounds.i2lin(iE);
            const Real Ye = YeBounds.i2lin(iYe);
            lambdas_d(iYe, iE, irho, 0) = Ye;
            diffs_t_d(iYe, iE, irho) =
                sc_d.TemperatureFromDensityInternalEnergy(
                    rho, sie, &(lambdas_d(iYe, iE, irho, 0)));
          });
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::fence();
#endif
      auto stop = std::chrono::high_resolution_clock::now();
      durationRoot += std::chrono::duration_cast<duration>(stop - start);
    }

    std::cout << "\t...Computing difference from ideal gas..." << std::endl;
    portableFor(
        "ideal gas T(rho, sie, Ye)", 0, nfine, 0, nfine, 0, nfine,
        PORTABLE_LAMBDA(const int iYe, const int iE, const int irho) {
          const Real rho = lRhoBounds.i2lin(irho);
          const Real sie = lEBounds.i2lin(iE);
          diffs_t_d(iYe, iE, irho) -=
              ig_d.TemperatureFromDensityInternalEnergy(rho, sie);
        });
#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::deep_copy(diffs_t_hv, diffs_t_dv);
#endif

    std::cout << "\t...Running once on host to histogram accesses..."
              << std::endl;
    std::vector<Real> lambdas_v(nfine*nfine*nfine*2);
    DataBox lambdas_h(lambdas_v.data(), nfine, nfine, nfine, 2);
    for (int trial = 0; trial < NTRIALS; ++trial) {
      for (int iYe = 0; iYe < nfine; ++iYe) {
        for (int iE = 0; iE < nfine; ++iE) {
          for (int irho = 0; irho < nfine; ++irho) {
            const Real rho = lRhoBounds.i2lin(irho);
            const Real sie = lEBounds.i2lin(iE);
            const Real Ye = YeBounds.i2lin(iYe);
            Real *lambda = &lambdas_h(iYe, iE, irho, 0);
            lambdas_h(iYe, iE, irho, 0) = Ye;
            sc.TemperatureFromDensityInternalEnergy(rho, sie, lambda);
          }
        }
      }
    }

    std::cout
        << "\nRESULTS:\n"
        << "Don't worry if errors are large. These are absolute, and include "
        << "extrapolation.\n"
        << "\tP(rho, T, Ye) time/point (microseconds)   = "
        << durationInterp.count() /
               static_cast<Real>(nfine * nfine * nfine * NTRIALS)
        << "\n"
        << "\tT(rho, sie, Ye) time/point (microseconds) = "
        << durationRoot.count() /
               static_cast<Real>(nfine * nfine * nfine * NTRIALS)
        << "\n"
        << "\tDelta P bounded by                        = "
        << std::max(std::abs(diffs_press_h.min()),
                    std::abs(diffs_press_h.max()))
        << "\n"
        << "\tDelta T bounded by                        = "
        << std::max(std::abs(diffs_t_h.min()), std::abs(diffs_t_h.max()))
        << "\n"
        << "Root finding:\n"
        << "its\tpercent taken:\n";
    Real tot = sc.counts.total();
    for (int i = 0; i < sc.counts.nBins(); ++i) {
      std::cout << i << "\t" << 100.0*sc.counts[i]/tot << "\n";
    }
    std::cout << std::endl;
    sc_d.Finalize();
    ig_d.Finalize();
  }
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::finalize();
#endif
  return 0;
}

#endif // SPINER_USE_HDF

//------------------------------------------------------------------------------
// © 2025. Triad National Security, LLC. All rights reserved.  This
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

// This example takes two materials and walks through valid sets of
// volume fractions alpha1, alpha2, such that the volume-averaged
// densities rhobar1 = alpha1 rho1 and rhobar2 = alpha2 rho2 are
// conserved and the volume fractions sum to 1. Here rho1 and rho2 are
// the microphysical densities of each material. At each volume
// fraction pair, it computes the temperature so that the total energy
// sums to a user specified value. When the pressures of each material
// agree at a given pair of volume fractions, pressure temperature
// equilibrium has been achived. This is useful for building an
// intuition for how pressure temperature equilibrium works and also
// seeing what it looks like for a given pairs of materials.  The
// script also outputs the residual as well as the P/T of each
// material at each point. This produces a "map" of the PT landscape
// for the material pair.

// C headers
#include <cmath>
#include <cstdio>

// C++ headers
#include <iostream>
#include <memory>
#include <string>

// This library contains portable utilities
#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_errors.hpp>

// This contains useful tools for preventing things like divide by zero
#include <singularity-eos/base/robust_utils.hpp>
// 1D root finding
#include <singularity-eos/base/root-finding-1d/root_finding.hpp>
// Needed to import the eos models
#include <singularity-eos/eos/eos.hpp>

// This library contains the spiner table object, which we will use to
// store our output
#include <spiner/databox.hpp>
#include <spiner/interpolation.hpp>

// These are the specializations of spiner we will use
using DataBox = Spiner::DataBox<Real>;
using RegularGrid1D = Spiner::RegularGrid1D<Real>;
using Spiner::DBDeleter;

// Number of equations of state to use, always 2
constexpr std::size_t NEOS = 2;

// Set the EOS you want to use here.
// Set the EOSs you want to use here.
using EOS = singularity::Variant<singularity::SpinerEOSDependsRhoT>;

int main(int argc, char *argv[]) {
  if (argc < 6) {
    std::cerr << "Usage: " << argv[0]
              << " residual_savename rhobar1 rhobar2 sietot nvfrac eosargs..."
              << std::endl;
    std::exit(1);
  }

  // This is needed for Kokkos
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::initialize();
#endif
  {
    // File we'll save the residual map to, in ascii. P/T trajectories will be printed to
    // to stdout
    const std::string savename = argv[1];
    // alpha rho for each material where here alpha is volume
    // fraction and rho is microphysical material density.
    const Real rhobar1 = atof(argv[2]);
    const Real rhobar2 = atof(argv[3]);
    // total specific internal energy in the problem.
    const Real sietot = atof(argv[4]);
    // number of vfrac points
    const std::size_t nvfrac = atoi(argv[5]);

    // Below are the arguments for an individual equation of
    // state. Here I'm assuming SpinerEOSDependsRhoT but you could
    // modify for your needs.
    const std::string materialsfile = argv[6];
    const int matids[] = {atoi(argv[7]), atoi(argv[8])};

    // Now let's load up the EOS's.
    EOS eos1 = singularity::SpinerEOSDependsRhoT(materialsfile, matids[0]);
    EOS eos2 = singularity::SpinerEOSDependsRhoT(materialsfile, matids[1]);

    // The total bulk density
    const Real rhotot = rhobar1 + rhobar2;
    const Real utot = rhotot * sietot;

    // We will be walking through volume fractions so let's grid that
    // up using a Spiner RegularGrid1D. This is the grid of alpha1s. alpha2 = 1 - alpha1.
    constexpr Real lvfrac_min = -5;
    constexpr Real lvfrac_max = -0.5;
    RegularGrid1D lvfracs(lvfrac_min, lvfrac_max, nvfrac);

    // We will be doing root finds in T space so we need to know what
    // our bounds are. There is often a minimum temperature for an EOS
    // but not a maximum.
    const Real Tmin = std::max(eos1.MinimumTemperature(), eos2.MinimumTemperature());
    const Real Tmax = 1e10; // change this if needed

    // We also want to store quantities of interest as we walk the volume fractions
    // This is a way of managing databoxes so that memory is automatically managed
    // Temperature
    std::unique_ptr<DataBox, DBDeleter> pTs(
        new DataBox(Spiner::AllocationTarget::Device, nvfrac));
    DataBox Ts = *pTs; // A shallow coppy of pTs. Unmanaged memory.
    // Pressure for materials 1 and 2
    std::unique_ptr<DataBox, DBDeleter> pPs(
        new DataBox(Spiner::AllocationTarget::Device, nvfrac, NEOS));
    DataBox Ps = *pPs;
    // vfrac, energy, and pressure residuals for the trajectory
    std::unique_ptr<DataBox, DBDeleter> pTrajRes(
        new DataBox(Spiner::AllocationTarget::Device, nvfrac, 3));
    DataBox trajRes = *pTrajRes;

    // We also need to ensure whatever internal tabulated data for
    // each EOS is available on device.
    eos1 = eos1.GetOnDevice();
    eos2 = eos2.GetOnDevice();

    // Ok let's walk the grid. We will do so on device with a portableFor loop
    portableFor(
        "Compute P and T for each alpha", 0, nvfrac, PORTABLE_LAMBDA(const int i) {
          // volume fractions
          const Real lvfrac = lvfracs.x(i);
          const Real vfrac = std::pow(10., lvfrac);
          const Real alpha1 = vfrac;
          const Real alpha2 = 1. - vfrac;

          // microphysical densities
          Real rho1 = singularity::robust::ratio(rhobar1, alpha1);
          Real rho2 = singularity::robust::ratio(rhobar2, alpha2);
          const Real Tguess = 1500; // just an initial guess

          // solve for the temperature such that sum of internal
          // energies is the total energy given microphysical densities rho
          auto fu = [=](const Real T) {
            return rhobar1 * eos1.InternalEnergyFromDensityTemperature(rho1, T) +
                   rhobar2 * eos2.InternalEnergyFromDensityTemperature(rho2, T);
          };
          // sets the T variable by reference
          auto root_status = RootFinding1D::regula_falsi(fu, utot, Tguess, Tmin, Tmax,
                                                         1e-16, 1e-16, Ts(i));
          if (root_status != RootFinding1D::Status::SUCCESS) {
            printf("# Root finder failed! "
                   "alpha1 alpha2, rho1, rho2, sie1(rho1,Tguess), sie2(rho2,Tguess) "
                   "sietot = "
                   "%.14e %.14e %.14e %.14e %.14e %.14e %.14e\n",
                   alpha1, alpha2, rho1, rho2,
                   eos1.InternalEnergyFromDensityTemperature(rho1, Tguess),
                   eos2.InternalEnergyFromDensityTemperature(rho2, Tguess), sietot);
            Ts(i) = Ps(i, 0) = Ps(i, 1) = -1;
            trajRes(i, 0) = trajRes(i, 1) = trajRes(i, 2) = -1;
          } else {
            // compute the pressures
            Ps(i, 0) = eos1.PressureFromDensityTemperature(rho1, Ts(i));
            Ps(i, 1) = eos2.PressureFromDensityTemperature(rho2, Ts(i));

            trajRes(i, 0) = alpha1 + alpha2 - 1.0;

            const Real u1 =
                rhobar1 * eos1.InternalEnergyFromDensityTemperature(rho1, Ts(i));
            const Real u2 =
                rhobar2 * eos2.InternalEnergyFromDensityTemperature(rho2, Ts(i));
            trajRes(i, 1) = singularity::robust::ratio(u1 + u2 - utot, utot);
            trajRes(i, 2) = singularity::robust::ratio(Ps(i, 1) - Ps(i, 0), utot);
          }
        });
    // wait for completion
#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::fence();
#endif

    // Copy the Ts array to host
    std::unique_ptr<DataBox, DBDeleter> pTs_h(
        new DataBox(Spiner::AllocationTarget::Host, nvfrac));
    DataBox Ts_h = *pTs_h;
    portableCopyToHost(Ts_h.data(), Ts.data(), Ts.sizeBytes());

    // Get max temperature achieved. Note this could have been done in
    // a single step via a portableReduce, but we split it out for clarity.
    Real Tmax_walked = Ts_h.max();

    // Let's make a T grid
    const Real nT = nvfrac + 1;

    RegularGrid1D gTs(Tmin, Tmax_walked, nT);

    // Now we can do one more loop where we compute the energy and
    // pressure residuals as a function of T and alpha
    std::unique_ptr<DataBox, DBDeleter> p_residuals(
        new DataBox(Spiner::AllocationTarget::Device, nT, nvfrac, NEOS));
    DataBox residuals = *p_residuals;

    portableFor(
        "Compute residual", 0, nT, 0, nvfrac, PORTABLE_LAMBDA(const int j, const int i) {
          // Temperature
          const Real T = gTs.x(j);

          // volume fractions
          const Real lvfrac = lvfracs.x(i);
          const Real vfrac = std::pow(10., lvfrac);
          const Real alpha1 = vfrac;
          const Real alpha2 = 1. - vfrac;

          // microphysical densities
          Real rho1 = singularity::robust::ratio(rhobar1, alpha1);
          Real rho2 = singularity::robust::ratio(rhobar2, alpha2);

          const Real u1 = rhobar1 * eos1.InternalEnergyFromDensityTemperature(rho1, T);
          const Real u2 = rhobar2 * eos2.InternalEnergyFromDensityTemperature(rho2, T);
          residuals(j, i, 0) = u1 + u2 - utot, utot;

          const Real P1 = eos1.PressureFromDensityTemperature(rho1, T);
          const Real P2 = eos2.PressureFromDensityTemperature(rho2, T);
          residuals(j, i, 1) = P1 - P2;
        });
#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::fence();
#endif

    // Now we need to copy all of this data to host so we can output it
    std::unique_ptr<DataBox, DBDeleter> pPs_h(
        new DataBox(Spiner::AllocationTarget::Host, nvfrac, NEOS));
    auto Ps_h = *pPs_h;
    portableCopyToHost(Ps_h.data(), Ps.data(), Ps.sizeBytes());

    std::unique_ptr<DataBox, DBDeleter> pTrajRes_h(
        new DataBox(Spiner::AllocationTarget::Host, nvfrac, 3));
    auto trajRes_h = *pTrajRes_h;
    portableCopyToHost(trajRes_h.data(), trajRes.data(), trajRes.sizeBytes());

    std::unique_ptr<DataBox, DBDeleter> p_residuals_h(
        new DataBox(Spiner::AllocationTarget::Host, nT, nvfrac, NEOS));
    auto residuals_h = *p_residuals_h;
    portableCopyToHost(residuals_h.data(), residuals.data(), residuals.sizeBytes());
#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::fence();
#endif

    FILE *file = fopen(savename.c_str(), "w");
    PORTABLE_REQUIRE(file != nullptr, "Unable to open file");
    // Now lets save the residuals
    fprintf(file, "# j i T alpha_1 res_sie res_P\n");
    for (std::size_t j = 0; j < nT; ++j) {
      for (std::size_t i = 0; i < nvfrac; ++i) {
        Real lvfrac = lvfracs.x(i);
        Real vfrac = std::pow(10., lvfrac);
        fprintf(file, "%ld %ld %.14e %.14e %.14e %.14e\n", j, i, gTs.x(j), vfrac,
                residuals_h(j, i, 0), residuals_h(j, i, 1));
      }
    }
    fclose(file);

    printf("# sie and P residuals normalized by total energy\n");
    printf("# i alpha_1 T P_1 P_2 res_alpha res_sie res_P\n");
    for (std::size_t i = 0; i < nvfrac; ++i) {
      Real lvfrac = lvfracs.x(i);
      Real vfrac = std::pow(10., lvfrac);
      printf("%ld %.14e %.14e %.14e %.14e %.14e %.14e %.14e\n", i, vfrac, Ts_h(i),
             Ps_h(i, 0), Ps_h(i, 1), trajRes(i, 0), trajRes(i, 1), trajRes(i, 2));
    }

    eos1.Finalize();
    eos2.Finalize();
  }
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::finalize();
#endif
  return 0;
}

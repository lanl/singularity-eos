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

#include <chrono>
#include <iostream>
#include <memory>
#include <stdlib.h>
#include <time.h>
#include <vector>

#include <ports-of-call/array.hpp>
#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>
#include <singularity-eos/closure/mixed_cell_models.hpp>
#include <spiner/databox.hpp>

#include <singularity-eos/eos/eos_models.hpp>
#include <singularity-eos/eos/eos_variant.hpp>

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#endif

#include <test/eos_unit_test_helpers.hpp>

using DataBox = Spiner::DataBox<Real>;
using singularity::ApproxTemperatureFromRhoMatU;
using singularity::IdealGas;
using singularity::PTESolverPT;
using singularity::PTESolverPTRequiredScratch;
using singularity::PTESolverRhoT;
using singularity::PTESolverRhoTRequiredScratch;
using singularity::Variant;
using EOS = Variant<IdealGas>;

template <typename RealIndexer, typename EOSIndexer>
PORTABLE_INLINE_FUNCTION Real set_state(Real rho_nom, Real sie_nom, RealIndexer &&rho,
                                        RealIndexer &&vfrac, RealIndexer &&sie,
                                        RealIndexer &&temp, RealIndexer &&press,
                                        EOSIndexer &&eos, std::size_t NMAT) {
  Real u_tot = 0;
  Real rho_tot = 0;
  for (std::size_t m = 0; m < NMAT; ++m) {
    // actual state
    rho[m] = rho_nom;
    sie[m] = sie_nom;

    // guesses
    temp[m] = 1;
    vfrac[m] = 1.0 / static_cast<Real>(NMAT);
    Real rhobar = vfrac[m] * rho[m];
    u_tot += rhobar * sie[m];
    rho_tot += rhobar;

    // pressure set
    // no lambdas... (CARE!)
    press[m] = eos[m].PressureFromDensityInternalEnergy(rho[m], sie[m]);
  }
  Real sie_tot = u_tot / rho_tot;
  return sie_tot;
}

SCENARIO("PT space PTE solver for two ideal gases", "[PTESolverPT][IdealGas]") {
  GIVEN("Two equations of state with different gammas") {
    constexpr std::size_t NEOS = 5;
    using EOS_Indexer_t = PortsOfCall::array<EOS, NEOS>;
    EOS_Indexer_t eoss;

    for (int m = 0; m < NEOS; ++m) {
      constexpr Real denom = NEOS + 1;
      eoss[m] = IdealGas((m + 1) / denom, m + 1);
    }

    WHEN("We create a thermodynamic state") {
      constexpr Real EPS = 1e-5;
      constexpr std::size_t NTRIAL = 1;
      const std::size_t nscratch_rt = PTESolverRhoTRequiredScratch(NEOS);
      const std::size_t nscratch_pt = PTESolverPTRequiredScratch(NEOS);
      const Real rho_nom = 1.0;
      const Real sie_nom = 1.0;
      int nsuccess = 0;
      portableReduce(
          "Test PTE ideal", 0, NTRIAL,
          PORTABLE_LAMBDA(const int t, int &s) {
            // Don't do this in a real problem!
            Real *rhos_save = (Real *)malloc(sizeof(Real) * NEOS);
            Real *rhos = (Real *)malloc(sizeof(Real) * NEOS);
            Real *vfracs = (Real *)malloc(sizeof(Real) * NEOS);
            Real *sies = (Real *)malloc(sizeof(Real) * NEOS);
            Real *Ts = (Real *)malloc(sizeof(Real) * NEOS);
            Real *Ps = (Real *)malloc(sizeof(Real) * NEOS);
            Real *lambda[NEOS];
            for (int i = 0; i < NEOS; i++) {
              lambda[i] = nullptr;
            }
            Real *scratch_rt = (Real *)malloc(sizeof(Real) * nscratch_rt);
            Real *scratch_pt = (Real *)malloc(sizeof(Real) * nscratch_pt);

            Real sie_tot =
                set_state(rho_nom, sie_nom, rhos, vfracs, sies, Ts, Ps, eoss, NEOS);
            Real rho_tot = 0;
            for (std::size_t m = 0; m < NEOS; ++m) {
              rho_tot += rhos[m] * vfracs[m];
            }
            const Real Tguess =
                ApproxTemperatureFromRhoMatU(NEOS, eoss, rho_tot * sie_tot, rhos, vfracs);

            auto method_rt = PTESolverRhoT<EOS_Indexer_t, Real *, decltype(lambda)>(
                NEOS, eoss, 1.0, sie_tot, rhos, vfracs, sies, Ts, Ps, lambda, scratch_rt,
                Tguess);
            auto status_rt = PTESolver(method_rt);
            bool success = status_rt.converged;

            for (std::size_t m = 0; m < NEOS; ++m) {
              rhos_save[m] = rhos[m];
            }

            auto method_pt = PTESolverPT<EOS_Indexer_t, Real *, decltype(lambda)>(
                NEOS, eoss, 1.0, sie_tot, rhos, vfracs, sies, Ts, Ps, lambda, scratch_pt,
                Tguess);
            auto status_pt = PTESolver(method_pt);
            success = success && status_pt.converged;

            for (std::size_t m = 0; m < NEOS; ++m) {
              bool matches = (std::abs(rhos[m] - rhos_save[m]) <= EPS);
              success = success && matches;
              if (!matches) {
                printf("rho doesn't match!\n"
                       "%ld %.14e %.14e %.14e\n",
                       m, rhos[m], rhos_save[m], std::abs(rhos[m] - rhos_save[m]));
              }
            }

            free(rhos_save);
            free(rhos);
            free(vfracs);
            free(sies);
            free(Ts);
            free(Ps);
            free(scratch_rt);
            free(scratch_pt);

            s += success;
          },
          nsuccess);

      THEN("A PTE solution can be found") { REQUIRE(nsuccess == NTRIAL); }
    }
  }
}

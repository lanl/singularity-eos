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

#include <chrono>
#include <iostream>
#include <memory>
#include <stdlib.h>
#include <time.h>
#include <vector>

// #include <ports-of-call/array.hpp>
#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>
#include <singularity-eos/base/indexable_types.hpp>

#include <singularity-eos/closure/mixed_cell_models.hpp>

#include <singularity-eos/eos/eos_models.hpp>
#include <singularity-eos/eos/eos_variant.hpp>

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#endif

#include <test/eos_unit_test_helpers.hpp>

using singularity::ApproxTemperatureFromRhoMatU;
using singularity::IdealElectrons;
using singularity::IdealGas;
using singularity::MeanAtomicProperties;
using singularity::PTESolverPT;
using singularity::PTESolverPTRequiredScratch;
using singularity::PTESolverRhoT;
using singularity::PTESolverRhoTRequiredScratch;
using singularity::Variant;
using EOS = Variant<IdealGas, IdealElectrons>;

using singularity::IndexableTypes::MeanIonizationState;
struct LambdaIndexerSingle {
  PORTABLE_FORCEINLINE_FUNCTION
  Real &operator[](const int i) { return z; }
  PORTABLE_FORCEINLINE_FUNCTION
  Real &operator[](const MeanIonizationState &s) { return z; }
  Real z = 0.9;
};

struct LambdaIndexer {
  PORTABLE_FORCEINLINE_FUNCTION
  auto &operator[](const int i) { return indexer; }
  PORTABLE_FORCEINLINE_FUNCTION
  const auto &operator[](const int i) const { return indexer; }
  LambdaIndexerSingle indexer;
};

template <typename RealIndexer, typename EOSIndexer>
PORTABLE_INLINE_FUNCTION Real set_state(Real rho_nom, Real sie_nom, RealIndexer &&rho,
                                        RealIndexer &&vfrac, RealIndexer &&sie,
                                        RealIndexer &&temp, RealIndexer &&press,
                                        EOSIndexer &&eos, std::size_t NMAT) {
  LambdaIndexer lambda;
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
    press[m] = eos[m].PressureFromDensityInternalEnergy(rho[m], sie[m], lambda);
  }
  Real sie_tot = u_tot / rho_tot;
  return sie_tot;
}

template <typename EOS_Indexer_t>
int ComparePTEs(EOS_Indexer_t eoss, const std::size_t NEOS, const std::size_t NTRIAL) {
  constexpr Real EPS = 1e-5;
  const std::size_t nscratch_rt = PTESolverRhoTRequiredScratch(NEOS, false);
  const std::size_t nscratch_pt = PTESolverPTRequiredScratch(NEOS, false);
  const Real rho_nom = 1.0;
  const Real sie_nom = 1.0;

  singularity::MixParams params;
  params.pte_rel_tolerance_e = 1e-14;
  params.pte_abs_tolerance_e = 1e-14;
  params.pte_abs_tolerance_v = 1e-14;
  params.pte_rel_tolerance_v = 1e-14;
  params.pte_rel_tolerance_p = 1e-14;
  params.pte_abs_tolerance_p = 1e-14;
  params.pte_residual_tolerance = 1.e-14;

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
        LambdaIndexer lambda;

        Real *scratch_rt = (Real *)malloc(sizeof(Real) * nscratch_rt);
        Real *scratch_pt = (Real *)malloc(sizeof(Real) * nscratch_pt);

        Real sie_tot =
            set_state(rho_nom, sie_nom, rhos, vfracs, sies, Ts, Ps, eoss, NEOS);
        Real rho_tot = 0;
        for (std::size_t m = 0; m < NEOS; ++m) {
          rho_tot += rhos[m] * vfracs[m];
        }
        const Real Tguess = ApproxTemperatureFromRhoMatU(NEOS, eoss, rho_tot * sie_tot,
                                                         rhos, vfracs, 0.0, lambda);

        auto method_rt = PTESolverRhoT<EOS_Indexer_t, Real *, decltype(lambda)>(
            NEOS, eoss, 1.0, sie_tot, rhos, vfracs, sies, Ts, Ps, lambda, scratch_rt,
            Tguess, params);
        auto status_rt = PTESolver(method_rt);
        bool success = status_rt.converged;

        for (std::size_t m = 0; m < NEOS; ++m) {
          rhos_save[m] = rhos[m];
        }

        sie_tot = set_state(rho_nom, sie_nom, rhos, vfracs, sies, Ts, Ps, eoss, NEOS);
        rho_tot = 0;
        for (std::size_t m = 0; m < NEOS; ++m) {
          rho_tot += rhos[m] * vfracs[m];
        }

        auto method_pt = PTESolverPT<EOS_Indexer_t, Real *, decltype(lambda)>(
            NEOS, eoss, 1.0, sie_tot, rhos, vfracs, sies, Ts, Ps, lambda, scratch_pt,
            Tguess, params);
        auto status_pt = PTESolver(method_pt);
        success = success && status_pt.converged;
        printf("PTE solvers converged in %ld iterations for RT and %ld iterations "
               "for PT\n",
               status_rt.max_niter, status_pt.max_niter);

        for (std::size_t m = 0; m < NEOS; ++m) {
          bool rho_matches = (std::abs(rhos[m] - rhos_save[m]) <= EPS);
          success = success && rho_matches;
          if (!rho_matches) {
            printf("rho doesn't match!\n"
                   "%ld %.14e %.14e %.14e\n",
                   m, rhos[m], rhos_save[m], std::abs(rhos[m] - rhos_save[m]));
          }
        }
        for (std::size_t m = 0; m < NEOS - 1; ++m) {
          bool P_matches = std::abs(Ps[m] - Ps[m + 1]) <=
                           0.5 * (std::abs(Ps[m]) + std::abs(Ps[m + 1]) + EPS) * EPS;
          if (!P_matches) {
            printf("Pressure not in equilibrium %ld %.14e %.14e!\n", m, Ps[m], Ps[m + 1]);
          }
          bool T_matches = std::abs(Ts[m] - Ts[m + 1]) <=
                           0.5 * (std::abs(Ts[m]) + std::abs(Ts[m + 1]) + EPS) * EPS;
          if (!T_matches) {
            printf("Temperatures not in equilibrium %ld %.14e %.14e!\n", m, Ts[m],
                   Ts[m + 1]);
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

  return nsuccess;
}

SCENARIO("PT space PTE solver for two ideal gases", "[PTESolverPT][IdealGas]") {
  constexpr std::size_t NEOS = 5;
  constexpr std::size_t NTRIAL = 1;
  using EOS_Indexer_t = std::array<EOS, NEOS>;
  GIVEN("Ideal gas equations of state with different gammas") {
    EOS_Indexer_t eoss;
    for (std::size_t m = 0; m < NEOS; ++m) {
      constexpr Real denom = NEOS + 1;
      eoss[m] = IdealGas((m + 1) / denom, m + 1);
    }
    WHEN("We compute PTE") {
      int nsuccess = ComparePTEs(eoss, NEOS, NTRIAL);
      THEN("A solution can be found") { REQUIRE(nsuccess == NTRIAL); }
    }
  }
  GIVEN("Electron equations of state with different mean atomic properties") {
    EOS_Indexer_t eoss;
    for (std::size_t m = 0; m < NEOS; ++m) {
      constexpr Real Abar = NEOS + 1;
      MeanAtomicProperties AZBar(Abar, m + 1);
      eoss[m] = IdealElectrons(AZBar);
    }
    WHEN("We compute PTE") {
      int nsuccess = ComparePTEs(eoss, NEOS, NTRIAL);
      THEN("A solution can be found") { REQUIRE(nsuccess == NTRIAL); }
    }
  }
}

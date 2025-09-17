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

#ifdef SINGULARITY_TEST_SESAME
#ifdef SINGULARITY_USE_SPINER_WITH_HDF5

#include <array>
#include <cmath>
#include <tuple>
#include <type_traits>
#include <utility>

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#endif

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_errors.hpp>

#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/indexable_types.hpp>
#include <singularity-eos/base/robust_utils.hpp>
#include <singularity-eos/base/variadic_utils.hpp>
#include <singularity-eos/eos/eos.hpp>
#include <singularity-eos/eos/eos_multi_eos.hpp>
#include <singularity-eos/eos/eos_spiner_rho_temp.hpp>
#include <test/eos_unit_test_helpers.hpp>

using singularity::make_MultiEOS;
using singularity::MassFracAverageFunctor;
using singularity::SpinerEOSDependsRhoT;
using singularity::Variant;
using singularity::VolumeFracHarmonicAverageFunctor;
using singularity::IndexableTypes::MassFraction;
using singularity::IndexerUtils::VariadicIndexer;

using Catch::Matchers::ContainsSubstring;

SCENARIO("Test the MultiEOS object with a binary alloy", "[MultiEOS][SpinerEOS]") {
  using namespace singularity;

  GIVEN("A pair of metal EOS for the binary alloy") {
    constexpr int num_eos = 2;
    constexpr int Al_matid = 3720;
    constexpr int Cu_matid = 3337;
    const std::string eos_file = "../materials.sp5";
    auto Cu_eos = SpinerEOSDependsRhoT{eos_file, Cu_matid};
    auto Al_eos = SpinerEOSDependsRhoT{eos_file, Al_matid};
    using LambdaT = VariadicIndexer<MassFraction<0>, MassFraction<1>>;

    // Lookup result tolerances
    constexpr Real lookup_tol = 1.0e-12;
    constexpr Real deriv_tol = lookup_tol * 1e3;
    constexpr Real dx_factor = 1.0e-08;
    constexpr Real finite_diff_tol = 1.0e-04;

    WHEN("The two EOS are combined in a MultiEOS object") {
      using BAvgF = VolumeFracHarmonicAverageFunctor;
      using GAvgF = VolumeFracHarmonicAverageFunctor;
      auto alloy = make_MultiEOS<BAvgF, GAvgF>(Cu_eos, Al_eos);

      // Set equal mass fractions
      std::array<Real, num_eos> set_mass_fracs{};
      set_mass_fracs.fill(1.0 / num_eos);
      LambdaT lambda{set_mass_fracs};
      const LambdaT lambda_init{set_mass_fracs};

      AND_WHEN("The alloy reference state values are returned") {
        Real rho_ref;
        Real temp_ref;
        Real sie_ref;
        Real pres_ref;
        Real cv_ref;
        Real bmod_ref;
        Real dpde_ref;
        Real dvdt_ref;
        alloy.ValuesAtReferenceState(rho_ref, temp_ref, sie_ref, pres_ref, cv_ref,
                                     bmod_ref, dpde_ref, dvdt_ref, lambda);

        // Sanity check to make sure the EOS are not overwriting mass fraction
        // values
        for (size_t i = 0; i < num_eos; i++) {
          INFO("i: " << i);
          REQUIRE_THAT(lambda[i], Catch::Matchers::WithinRel(lambda_init[i], lookup_tol));
        }

        INFO("pres_ref = " << pres_ref << "   temp_ref = " << temp_ref);

        THEN("They agree with the states from reference pressure and temperature") {
          std::array<Real, num_eos> density_mat;
          std::array<Real, num_eos> sie_mat;
          Real rho_PT;
          Real sie_PT;
          alloy.GetStatesFromPressureTemperature(rho_PT, sie_PT, pres_ref, temp_ref,
                                                 density_mat, sie_mat, lambda);

          // Sanity check to make sure the EOS are not overwriting mass fraction
          // values
          for (size_t i = 0; i < num_eos; i++) {
            INFO("i: " << i);
            REQUIRE_THAT(lambda[i],
                         Catch::Matchers::WithinRel(lambda_init[i], lookup_tol));
          }

          CHECK_THAT(rho_ref, Catch::Matchers::WithinRel(rho_PT, lookup_tol));
          CHECK_THAT(sie_ref, Catch::Matchers::WithinRel(sie_PT, lookup_tol));

          AND_THEN("The lookups on the individual EOS are consistent") {
            const auto eos_arr = alloy.CreateEOSArray();
            for (size_t i = 0; i < num_eos; i++) {
              INFO("EOS[" << i << "]: " << eos_arr[i].EosType());
              Real rho_mat_PT;
              Real sie_mat_PT;
              eos_arr[i].DensityEnergyFromPressureTemperature(pres_ref, temp_ref, lambda,
                                                              rho_mat_PT, sie_mat_PT);

              // Sanity check to make sure the EOS are not overwriting mass fraction
              // values
              for (size_t i = 0; i < num_eos; i++) {
                INFO("i: " << i);
                REQUIRE_THAT(lambda[i],
                             Catch::Matchers::WithinRel(lambda_init[i], lookup_tol));
              }

              CHECK_THAT(rho_mat_PT,
                         Catch::Matchers::WithinRel(density_mat[i], lookup_tol));
              CHECK_THAT(sie_mat_PT, Catch::Matchers::WithinRel(sie_mat[i], lookup_tol));
            }
          }
        }

        THEN("They agree with values from FillEos at the reference pressure and "
             "temperature") {
          constexpr unsigned long input = thermalqs::pressure | thermalqs::temperature;
          constexpr unsigned long output = ~input;
          Real pres_FillEos = pres_ref;
          Real temp_FillEos = temp_ref;
          Real rho_FillEos;
          Real sie_FillEos;
          Real cv_FillEos;
          Real bmod_FillEos;
          alloy.FillEos(rho_FillEos, temp_FillEos, sie_FillEos, pres_FillEos, cv_FillEos,
                        bmod_FillEos, output, lambda);

          // Sanity check to make sure the EOS are not overwriting mass fraction
          // values
          for (size_t i = 0; i < num_eos; i++) {
            INFO("i: " << i);
            REQUIRE_THAT(lambda[i],
                         Catch::Matchers::WithinRel(lambda_init[i], lookup_tol));
          }

          CHECK_THAT(rho_ref, Catch::Matchers::WithinRel(rho_FillEos, lookup_tol));
          CHECK_THAT(sie_ref, Catch::Matchers::WithinRel(sie_FillEos, lookup_tol));
          CHECK_THAT(cv_ref, Catch::Matchers::WithinRel(cv_FillEos, deriv_tol));
          CHECK_THAT(bmod_ref, Catch::Matchers::WithinRel(bmod_FillEos, deriv_tol));
        }

        THEN("dpde and dvdt agree with finite difference approximations") {
          const Real dT = temp_ref * dx_factor + dx_factor;
          INFO("  dT = " << dT);

          Real rho_pert;
          Real sie_pert;
          alloy.DensityEnergyFromPressureTemperature(pres_ref, temp_ref + dT, lambda,
                                                     rho_pert, sie_pert);

          // Sanity check to make sure the EOS are not overwriting mass fraction
          // values
          for (size_t i = 0; i < num_eos; i++) {
            INFO("i: " << i);
            REQUIRE_THAT(lambda[i],
                         Catch::Matchers::WithinRel(lambda_init[i], lookup_tol));
          }

          const Real dV = robust::ratio(1.0, rho_pert) - robust::ratio(1.0, rho_ref);
          const Real dvdt_FD = dV / dT;

          INFO("  dV = " << dV);
          CHECK_THAT(dvdt_ref, Catch::Matchers::WithinRel(dvdt_FD, finite_diff_tol));

          // Calculate a de at constant volume from dT since density-energy PTE solver
          // can require looser tolerances
          sie_pert =
              alloy.InternalEnergyFromDensityTemperature(rho_ref, temp_ref + dT, lambda);
          const Real pres_pert =
              alloy.PressureFromDensityTemperature(rho_ref, temp_ref + dT, lambda);

          // Sanity check to make sure the EOS are not overwriting mass fraction
          // values
          for (size_t i = 0; i < num_eos; i++) {
            INFO("i: " << i);
            REQUIRE_THAT(lambda[i],
                         Catch::Matchers::WithinRel(lambda_init[i], lookup_tol));
          }

          const Real dP = pres_pert - pres_ref;
          const Real de = sie_pert - sie_ref;
          const Real dpde_FD = dP / de;
          INFO("  de (from dT) = " << de << "  dP = " << dP);

          CHECK_THAT(dpde_ref, Catch::Matchers::WithinRel(dpde_FD, finite_diff_tol));
        }
      }
    }
  }
}

#endif // SINGULARITY_USE_SPINER_WITH_HDF5
#endif // SINGULARITY_TEST_SESAME

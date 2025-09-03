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

#include <array>
#include <cmath>
#include <tuple>
#include <vector>

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#endif

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_errors.hpp>

#include <singularity-eos/base/indexable_types.hpp>
#include <singularity-eos/base/robust_utils.hpp>
#include <singularity-eos/eos/eos.hpp>
#include <singularity-eos/eos/eos_multi_eos.hpp>
#include <test/eos_unit_test_helpers.hpp>

using singularity::DavisProducts;
using singularity::DavisReactants;
using singularity::ShiftedEOS;
using singularity::Variant;
using singularity::IndexableTypes::MassFraction;
using singularity::IndexerUtils::VariadicIndexer;

using Catch::Matchers::ContainsSubstring;

SCENARIO("Test the MultiEOS object with reactants and products EOS",
         "[MultiEOS][DavisReactants][DavisProducts]") {
  GIVEN("A pair of EOS for reactants and products and a resulting MultiEOS "
        "object") {
    // Unit conversions
    [[maybe_unused]] constexpr Real cm = 1.;
    [[maybe_unused]] constexpr Real us = 1e-06;
    [[maybe_unused]] constexpr Real Mbcc_per_g = 1e12;
    constexpr Real GPa = 1.0e10;
    constexpr Real sqrtGPa = 1.0e5; // sqrt isn't constexpr...
    constexpr Real MJ_per_kg = 1.0e10;

    // Davis Reactants EOS
    constexpr Real rho0_DP = 1.890;   // g/cm^3
    constexpr Real e0_DP = 0.;        // erg / g
    constexpr Real P0_DP = 0.;        // microbar
    constexpr Real T0_DP = 297;       // K
    constexpr Real A = 1.8 * sqrtGPa; // sqrt(microbar)
    constexpr Real B = 4.6;
    constexpr Real C = 0.34;
    constexpr Real G0 = 0.56;
    constexpr Real Z = 0.0;
    constexpr Real alpha = 0.4265;
    constexpr Real Cv_DP = 0.001074 * MJ_per_kg; // erg / g / K
    auto davis_r_eos =
        DavisReactants(rho0_DP, e0_DP, P0_DP, T0_DP, A, B, C, G0, Z, alpha, Cv_DP);

    // Davis Products EOS
    constexpr Real a = 0.798311;
    constexpr Real b = 0.58;
    constexpr Real k = 1.35;
    constexpr Real n = 2.66182;
    constexpr Real vc = 0.75419;              // cm^3 / g
    constexpr Real pc = 3.2 * GPa;            // microbar
    constexpr Real Cv = 0.001072 * MJ_per_kg; // erg / g / K
    constexpr Real Q = 4.115 * MJ_per_kg;
    auto davis_p_eos = DavisProducts(a, b, k, n, vc, pc, Cv).Modify<ShiftedEOS>(-Q);

    // Create the multiEOS object
    constexpr size_t num_eos = 2;
    using LambdaT = VariadicIndexer<MassFraction<0>, MassFraction<1>>;
    auto multi_eos = make_MultiEOS(davis_r_eos, davis_p_eos);

    WHEN("A mass fraction cutoff less than 0 is specified") {
      [[maybe_unused]] constexpr Real mf_cutoff = -0.01;
      THEN("Constructing the MultiEOS object should throw an exception") {
        // REQUIRE_MAYBE_THROWS_WITH(make_MultiEOS(mf_cutoff, davis_r_eos, davis_p_eos),
        //                           "Mass fracton cutoff must be non-negative");
        REQUIRE_MAYBE_THROWS(make_MultiEOS(mf_cutoff, davis_r_eos, davis_p_eos));
      }
    }

    WHEN("A mass fraction cutoff equal to 1 is specified") {
      [[maybe_unused]] constexpr Real mf_cutoff = 1.0;
      THEN("Constructing the MultiEOS object should throw an exception") {
        // REQUIRE_MAYBE_THROWS_WITH(make_MultiEOS(mf_cutoff, davis_r_eos, davis_p_eos),
        //                           "Mass fractons must be less than 1");
        REQUIRE_MAYBE_THROWS(make_MultiEOS(mf_cutoff, davis_r_eos, davis_p_eos));
      }
    }

    WHEN("A mass fraction cutoff greater than 1 is specified") {
      [[maybe_unused]] constexpr Real mf_cutoff = 2.0;
      THEN("Constructing the MultiEOS object should throw an exception") {
        // REQUIRE_MAYBE_THROWS_WITH(make_MultiEOS(mf_cutoff, davis_r_eos, davis_p_eos),
        //                           "Mass fractons must be less than 1");
        REQUIRE_MAYBE_THROWS(make_MultiEOS(mf_cutoff, davis_r_eos, davis_p_eos));
      }
    }

    WHEN("A set of mass fractions is provided in a lambda") {
      std::array<Real, num_eos> set_mass_fracs{};
      set_mass_fracs.fill(1.0 / num_eos);
      LambdaT lambda{set_mass_fracs};

      THEN("The MultiEOS object can translate these from the lambda to an array") {
        std::array<Real, num_eos> mass_fracs{};
        multi_eos.PopulateMassFracArray(mass_fracs, lambda);
        for (size_t m = 0; m < num_eos; m++) {
          INFO("Mass fraction test: mass_fracs[" << m << "] = " << mass_fracs[m]);
          CHECK_THAT(set_mass_fracs[m],
                     Catch::Matchers::WithinRel(mass_fracs[m], 1.0e-14));
        }
      }
    }

    THEN("The MultiEOS object can be placed in an EOS Variant") {
      using EOS =
          Variant<DavisReactants, ShiftedEOS<DavisProducts>,
                  decltype(make_MultiEOS(DavisReactants{}, ShiftedEOS<DavisProducts>{}))>;
      EOS multi_eos_in_variant = make_MultiEOS(davis_r_eos, davis_p_eos);

      AND_WHEN("A pressure-temperature lookup is performed") {
        // Populate lambda with mass fractions (only)
        std::array<Real, num_eos> set_mass_fracs{};
        set_mass_fracs.fill(1.0 / num_eos);
        LambdaT lambda{set_mass_fracs};

        // A high pressure and temperature
        constexpr Real P = 8.0e10;
        constexpr Real T = 8000;
        Real rho;
        Real sie;
        multi_eos_in_variant.DensityEnergyFromPressureTemperature(P, T, lambda, rho, sie);
        INFO("MultiEOS density: " << rho << "   sie: " << sie);
        CHECK(rho > 0.);
        CHECK(std::isnormal(rho));
        CHECK(std::isnormal(sie));

        THEN("The result is consistent with doing the same for the individual EOS") {

          // Create an array of EOS
          const auto eos_arr = multi_eos.CreateEOSArray();

          // Create an array of mass fractions for adding up the individual values
          std::array<Real, num_eos> mass_fracs{};
          multi_eos.PopulateMassFracArray(mass_fracs, lambda);

          // Sum the individual EOS contributions weighted by mass fractions
          Real spvol_bulk = 0.;
          Real sie_bulk = 0.;
          for (size_t m = 0; m < num_eos; m++) {
            Real rho_m;
            Real sie_m;
            eos_arr[m].DensityEnergyFromPressureTemperature(P, T, nullptr, rho_m, sie_m);
            INFO("EOS[" << m << "]:\n" << eos_arr[m].EosType());
            CHECK(rho_m > 0.);
            CHECK(sie_m != 0.);
            spvol_bulk += mass_fracs[m] / rho_m;
            sie_bulk += mass_fracs[m] * sie_m;
          }

          // Quick check that we can invert the specific volume
          REQUIRE(spvol_bulk > 0.);
          Real rho_bulk = 1.0 / spvol_bulk;

          // Make sure the weighted sums add up to the result from the MultiEOS
          CHECK_THAT(rho_bulk, Catch::Matchers::WithinRel(rho, 1.0e-12));
          CHECK_THAT(sie_bulk, Catch::Matchers::WithinRel(sie, 1.0e-12));
        }
      }
    }

    WHEN("A mass fraction cutoff is specified in the constructor") {
      constexpr Real mf_cutoff = 0.001;
      auto multi_eos_largeMF = make_MultiEOS(mf_cutoff, davis_r_eos, davis_p_eos);

      THEN("The 'SetUpPTE' member function will correctly determine which materials are "
           "participating in PTE") {
        constexpr auto small_mf = mf_cutoff / 10.;
        constexpr auto large_mf = (1.0 - small_mf) / (num_eos - 1);

        // Somewhat realistic but irrelvant values
        constexpr Real density_tot = 1.5;
        constexpr Real sie_tot = 1.0e10;

        // Arrays
        std::array<Real, num_eos> mass_fracs{};
        mass_fracs.fill(large_mf);
        std::array<Real, num_eos> vol_fracs{};   // zero initialized
        std::array<Real, num_eos> density_mat{}; // zero initialized
        std::array<Real, num_eos> sie_mat{};     // zero initialized
        std::array<size_t, num_eos> pte_mats{};  // zero initialized
        size_t npte = 0;
        Real vfrac_tot = 0;
        auto eos_arr = multi_eos_largeMF.CreateEOSArray();

        // Select first index to have a small mass fraction so that it won't
        // participate
        mass_fracs[0] = small_mf;

        multi_eos_largeMF.SetUpPTE(mass_fracs, vol_fracs, density_mat, sie_mat, pte_mats,
                                   npte, vfrac_tot, density_tot, sie_tot, eos_arr);

        INFO("mass_fracs[0] = " << mass_fracs[0]);
        INFO("Small Mass Fraction = " << small_mf);
        INFO("Mass Fraction Cutoff = " << mf_cutoff);
        CHECK(npte == num_eos - 1);
        CHECK(pte_mats[0] != 0); // First PTE material should not be index 0
      }
    }

    WHEN("The EosType function is called") {
      auto davis_r_eostype = davis_r_eos.EosType();
      auto davis_p_eostype = davis_r_eos.EosType();
      auto multi_eostype = multi_eos.EosType();

      THEN("The MultiEOS EosType should contain the component EOS names") {
        CHECK_THAT(multi_eostype, ContainsSubstring(davis_r_eostype));
        CHECK_THAT(multi_eostype, ContainsSubstring(davis_p_eostype));
        CHECK_THAT(multi_eostype, ContainsSubstring("MultiEOS"));
      }
    }

    WHEN("The EosPyType function is called") {
      auto davis_r_eospytype = davis_r_eos.EosPyType();
      auto davis_p_eospytype = davis_r_eos.EosPyType();
      auto multi_eospytype = multi_eos.EosPyType();

      THEN("The MultiEOS EosType should contain the component EOS names") {
        CHECK_THAT(multi_eospytype, ContainsSubstring(davis_r_eospytype));
        CHECK_THAT(multi_eospytype, ContainsSubstring(davis_p_eospytype));
        CHECK_THAT(multi_eospytype, ContainsSubstring("MultiEOS"));
      }
    }
  }
}

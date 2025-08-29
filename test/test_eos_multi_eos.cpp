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
#include <catch2/matchers/catch_matchers_string.hpp>
#endif

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_errors.hpp>

#include <singularity-eos/base/robust_utils.hpp>
#include <singularity-eos/eos/eos.hpp>
#include <singularity-eos/eos/eos_multi_eos.hpp>
#include <test/eos_unit_test_helpers.hpp>

using singularity::DavisProducts;
using singularity::DavisReactants;
using singularity::ShiftedEOS;
using singularity::Variant;

using Catch::Matchers::ContainsSubstring;

SCENARIO("Test the MultiEOS object with reactants and products EOS",
         "[MultiEOS][DavisReactants][DavisProducts]") {
  GIVEN("A pair of EOS for reactants and products and a resulting MultiEOS "
        "object") {
    // Unit conversions
    constexpr Real cm = 1.;
    constexpr Real us = 1e-06;
    constexpr Real Mbcc_per_g = 1e12;
    constexpr Real GPa = 1.0e10;
    constexpr Real MJ_per_kg = 1.0e10;

    // Davis Reactants EOS
    constexpr Real rho0_DP = 1.890;          // g/cm^3
    constexpr Real e0_DP = 0.;               // erg / g
    constexpr Real P0_DP = 0.;               // microbar
    constexpr Real T0_DP = 297;              // K
    constexpr Real A = 1.8 * std::sqrt(GPa); // sqrt(microbar)
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
    auto multi_eos = make_MultiEOS(davis_r_eos, davis_p_eos);

    THEN("The MultiEOS object can also be placed in an EOS Variant") {
      using EOS =
          Variant<DavisReactants, ShiftedEOS<DavisProducts>,
                  decltype(make_MultiEOS(DavisReactants{}, ShiftedEOS<DavisProducts>{}))>;
      EOS multi_eos_in_variant = make_MultiEOS(davis_r_eos, davis_p_eos);
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

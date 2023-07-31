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

#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#ifndef CATCH_CONFIG_RUNNER
#include "catch2/catch.hpp"
#endif

#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/eos/eos.hpp>
#include <test/eos_unit_test_helpers.hpp>

using singularity::EOS;
using singularity::SAP_Polynomial;

SCENARIO("SAP_Polynomial EOS", "Check if eos returns expected values") {
  GIVEN("Parameters for a SAP_Polynomial EOS") {
    // Unit conversions
    constexpr Real ubar_per_Mbar = 1e12;

    // Constants of the form 0.0? are zero for the original eos
    constexpr Real rho0 = 16.6;
    constexpr Real a0 = 0.01 * ubar_per_Mbar;
    constexpr Real a1 = 1.09 * ubar_per_Mbar;
    constexpr Real a2c = 0.90 * ubar_per_Mbar;
    constexpr Real a2e = 0.02 * ubar_per_Mbar;
    constexpr Real a3 = 1.80 * ubar_per_Mbar;
    constexpr Real b0 = 1.90 * ubar_per_Mbar;
    constexpr Real b1 = 0.01 * ubar_per_Mbar;
    constexpr Real b2c = 0.02 * ubar_per_Mbar;
    constexpr Real b2e = 0.03 * ubar_per_Mbar;
    constexpr Real b3 = 0.04 * ubar_per_Mbar;

    // Create the EOS
    EOS host_eos = SAP_Polynomial(rho0, a0, a1, a2c, a2e, a3, b0, b1, b2c, b2e, b3);
    EOS eos = host_eos.GetOnDevice();

    GIVEN("Densities and energies") {

      // Size of input
      constexpr int num = 5;

      // Arrays for input/output
#ifdef PORTABILITY_STRATEGY_KOKKOS
      // Create Kokkos views on device. Output arrays will have mirror views
      Kokkos::View<Real[num]> v_rho("rho");
      Kokkos::View<Real[num]> v_sie("sie");

      Kokkos::View<Real[num]> v_P("Pressure");
      Kokkos::View<Real[num]> v_B("Bulk modulus");
      Kokkos::View<Real[num]> v_gamma("Gruneisen parameter");
      auto P = Kokkos::create_mirror_view(v_P);
      auto B = Kokkos::create_mirror_view(v_B);
      auto gamma = Kokkos::create_mirror_view(v_gamma);
#else
      // Otherwise create arrays to contain values and create pointers to
      // be passed to the functions in place of the Kokkos views
      std::array<Real, num> rho;
      std::array<Real, num> sie;
      std::array<Real, num> P;
      std::array<Real, num> B;
      std::array<Real, num> gamma;
      auto v_rho = rho.data();
      auto v_sie = sie.data();
      auto v_P = P.data();
      auto v_B = B.data();
      auto v_gamma = gamma.data();
#endif // PORTABILITY_STRATEGY_KOKKOS

      // Fill input arrays
      portableFor(
          "Initialize rho, siw", 0, 1, PORTABLE_LAMBDA(int i) {
            v_rho[0] = 15.4;
            v_rho[1] = 16.0;
            v_rho[2] = 16.6;
            v_rho[3] = 17.2;
            v_rho[4] = 17.8;

            v_sie[0] = 1e-3;
            v_sie[1] = 1e-3;
            v_sie[2] = 1e-3;
            v_sie[3] = 1e-3;
            v_sie[4] = 1e-3;
          });

      // Expected values from calculations
      constexpr std::array<Real, num> P_expected{
          -6.747122099662988281250000000000e+10, -2.755678257812800216674804687500e+10,
          1.190000000000000000000000000000e+10, 5.255876399778214263916015625000e+10,
          9.607914667524797058105468750000e+10};

      constexpr std::array<Real, num> B_expected{
          1.034707096191414062500000000000e+12, 1.056016317964556762695312500000e+12,
          1.090010000000000122070312500000e+12, 1.204131143205424072265625000000e+12,
          1.338595278608992675781250000000e+12};

      constexpr std::array<Real, num> gamma_expected{
          1.899418769576782958984375000000e+12, 1.899675858317870117187500000000e+12,
          1.900000000000000000000000000000e+12, 1.900389463209202148437500000000e+12,
          1.900842516531505615234375000000e+12};

#ifdef PORTABILITY_STRATEGY_KOKKOS
      // Create host-side mirrors of the inputs and copy the inputs. These are
      // just used for the comparisons
      auto rho = Kokkos::create_mirror_view(v_rho);
      auto sie = Kokkos::create_mirror_view(v_sie);
      Kokkos::deep_copy(rho, v_rho);
      Kokkos::deep_copy(sie, v_sie);
#endif // PORTABILITY_STRATEGY_KOKKOS

      WHEN("A P(rho, e) lookup is performed") {
        eos.PressureFromDensityInternalEnergy(v_rho, v_sie, v_P, num);
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(P, v_P);
#endif
        THEN("The returned P(rho, e) should be equal to the true value") {
          array_compare(num, rho, sie, P, P_expected, "Density", "Energy");
        }
      }

      WHEN("A B(rho, e) lookup is performed") {
        eos.BulkModulusFromDensityInternalEnergy(v_rho, v_sie, v_B, num);
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(B, v_B);
#endif
        THEN("The returned B(rho, e) should be equal to the true value") {
          array_compare(num, rho, sie, B, B_expected, "Density", "Energy");
        }
      }

      WHEN("A gamma(rho, e) lookup is performed") {
        eos.GruneisenParamFromDensityInternalEnergy(v_rho, v_sie, v_gamma, num);
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(gamma, v_gamma);
#endif
        THEN("The returned gamma(rho, e) should be equal to the true value") {
          array_compare(num, rho, sie, gamma, gamma_expected, "Density", "Energy");
        }
      }
    }
  }
}

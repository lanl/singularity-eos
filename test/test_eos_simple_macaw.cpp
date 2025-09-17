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

#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#endif

#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/eos/eos.hpp>
#include <test/eos_unit_test_helpers.hpp>

using singularity::SimpleMACAW;
using EOS = singularity::Variant<SimpleMACAW>;

constexpr Real REAL_TOL = 1e-12;

// Only run exception checking when we aren't offloading to the GPUs
SCENARIO("Testing the Simple MACAW EOS", "[SimpleMACAWEOS]") {
  GIVEN("Parameters for a Simple MACAW EOS") {
    // Parameters for copper (See Table 1 in the paper by Aslam & Lozano)
    constexpr Real A = 7.3;
    constexpr Real B = 3.9;
    constexpr Real Cvinf = 0.000389;
    constexpr Real v0 = 1. / 8.952;
    constexpr Real T0 = 150.;
    constexpr Real Gc = 0.5;

    // Create the EOS
    auto host_eos = SimpleMACAW(A, B, Cvinf, v0, T0, Gc);
    auto eos = host_eos.GetOnDevice();

    WHEN("A set of densities is provided") {
      Real rho = 0.5;
      for (int i = 0; i < 10; i++) {
        rho += rho + i; // cylce through a variety of densities
        DYNAMIC_SECTION("For a given density " << rho << " and energy from the cold curve") {
          Real e = eos.SieColdCurve(rho);
          INFO("rho = " << rho << "  e = " << e);
          THEN("The temperature at this density and energy should be zero") {
            REQUIRE(eos.TemperatureFromDensityInternalEnergy(rho, e) == 0.0);
          } // Then
        } // Dynamic Section
      } // for
    } // When
    host_eos.Finalize();
    eos.Finalize();
  } // Given
} // Scenario

SCENARIO("Testing the Variant API Simple MACAW EOS", "[SimpleMACAWEOS]") {
  GIVEN("Parameters for a Simple MACAW EOS") {
    // Parameters for copper (See Table 1 in the paper by Aslam & Lozano)
    constexpr Real A = 7.3;
    constexpr Real B = 3.9;
    constexpr Real Cvinf = 0.000389;
    constexpr Real v0 = 1. / 8.952;
    constexpr Real T0 = 150.;
    constexpr Real Gc = 0.5;

    // Create the EOS
    using EOS = singularity::Variant<SimpleMACAW>;
    EOS eos_in_variant = SimpleMACAW(A, B, Cvinf, v0, T0, Gc);

    WHEN("Densities are provided") {
      // For large temperatures, the specific heat capacity should approach a constant value
      // (The Dulong-Petit Law)
      Real rho = 0.1;
      Real T = 1e7;
      Real Cv = eos_in_variant.SpecificHeatFromDensityTemperature(rho, T);
      for (int i = 0; i < 10; i++) {
        rho += rho + i; // cylce through a variety of densities
        DYNAMIC_SECTION("For a given density " << rho << " and temperature " << T) {
          INFO(std::fixed << std::setprecision(15) << 
               "Cv = " << Cv << ", Cvinf = " << Cvinf << ",  T = " << T << ", rho = " << rho);
          THEN("The Dulong-Petit law should be satisfied") {
            REQUIRE(isClose(Cv, Cvinf, 1e-8 * Cvinf));
          } // Then
        } // Dynamic section
      } // for

      // Check for zero specific heat at zero temperature
      Cv = eos_in_variant.SpecificHeatFromDensityTemperature(rho, 0.);
      THEN("The specific heat capacity at constant volume should be zero at zero temperature") {
        INFO(std::fixed << std::setprecision(15) << "Cv = " << Cv << ", rho = " << rho);
        REQUIRE(Cv == 0.0);
      }
    }

    // Check that the EOS is inverting correctly
    WHEN("A set of densities and energies are provided") {
      using EOS = singularity::Variant<SimpleMACAW>;
      EOS eos_in_variant = SimpleMACAW(A, B, Cvinf, v0, T0, Gc);

      const Real rho_0 = 0.3;
      const Real sie_0 = 190.0;
      const Real P = eos_in_variant.PressureFromDensityInternalEnergy(rho_0, sie_0);
      const Real T = eos_in_variant.TemperatureFromDensityInternalEnergy(rho_0, sie_0);
      Real rho, sie;
      Real* lambda;
      eos_in_variant.DensityEnergyFromPressureTemperature(P, T, lambda, rho, sie);
      INFO("rho_0 = " << rho_0 << ", rho = " << rho << ", sie_0 = " << sie_0 << ", sie = " << sie);
      INFO("P = " << P << "  T = " << T);
      THEN("The inversion back to rho and sie from P and T should be identital.") {
        REQUIRE(isClose(rho, rho_0, REAL_TOL * rho_0));
        REQUIRE(isClose(sie, sie_0, 100 * REAL_TOL * sie_0));
      }
    } // When
    eos_in_variant.Finalize();
  } // Given
} // Scenario

//------------------------------------------------------------------------------
// © 2021-2024. Triad National Security, LLC. All rights reserved.  This
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

constexpr Real REAL_TOL = std::numeric_limits<Real>::epsilon() * 1e4; // 1e-12 for double



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
          Real v = 1.0 / rho;
          Real e = eos.SieColdCurve(v);
          INFO("rho = " << rho << "  v = " << v << "  e = " << e);
          THEN("The temperature at this density and energy should be zero") {
            REQUIRE(eos.TemperatureFromDensityInternalEnergy(rho, e) == 0.0);
          } // Then
        } // Dynamic Section
      } // for
    } // When

    WHEN("A set of densities and energies are provided") {
      const Real rho_0 = 0.3;
      const Real sie_0 = 2.7;
      const Real P = eos.PressureFromDensityInternalEnergy(rho_0, sie_0);
      const Real T = eos.TemperatureFromDensityInternalEnergy(rho_0, sie_0);
      Real rho, sie;
      Real* lambda;
      eos.DensityEnergyFromPressureTemperature(P, T, lambda, rho, sie);
      sie = eos.InternalEnergyFromDensityPressure(rho, P);
      INFO("rho_0 = " << rho_0 << "  rho = " << rho << "  e_0 = " << sie_0 << "  e = " << sie);
      INFO("P = " << P << "  T = " << T);
      THEN("The inversion back to rho and sie from P and T should be identital.") {
        REQUIRE(isClose(rho, rho_0, REAL_TOL * rho_0));
        REQUIRE(isClose(sie, sie_0, REAL_TOL * sie_0));
      }
    } // When

    WHEN("") {
      const Real rho_0 = 3.;
      const Real T_0 = 178.;
      const Real Bs = eos.BulkModulusFromDensityTemperature(rho_0, T_0);
      const Real BT = eos.IsothermalBulkModulusFromDensityTemperature(rho_0, T_0);
      const Real cv = eos.SpecificHeatFromDensityTemperature(rho_0, T_0);
      const Real cp = eos.ConstantPressureSpecificHeatFromDensityTemperature(rho_0, T_0);
      THEN("Thermodynamic stability should satisfy: Bs >= BT and cp >= cv.") {
        REQUIRE(Bs >= BT);
        REQUIRE(cp >= cv);
      }
    }
  } // Given
} // Scenario

SCENARIO("Thermodynamic consistency", "[SimpleMACAWEOS][Thermo Consistency]") {
  GIVEN("Initial density and temperature") { 
  }
}

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

// Parameters for copper
constexpr Real A = 7.3;
constexpr Real B = 3.9;
constexpr Real Cvinf = 0.000389;
constexpr Real v0 = 1. / 8.952;
constexpr Real T0 = 150.;
constexpr Real Gc = 0.5;

// Create the EOS
auto host_eos = SimpleMACAW(A, B, Cvinf, v0, T0, Gc);
auto eos = host_eos.GetOnDevice();

// Only run exception checking when we aren't offloading to the GPUs
SCENARIO("Zero temperature on cold curve", "[SimpleMACAWEOS][Temperature]") {
  GIVEN("Parameters for a Simple MACAW EOS") {
    Real rho = 0.5;
    for (int i = 0; i < 10; i++) {
      rho += rho + i; // cylce through a variety of densities
      THEN("The temperature evaluated on the cold curve should produce zero") {
        Real v = 1.0 / rho;
        Real e = eos.SieColdCurve(v);
        REQUIRE(eos.TemperatureFromDensityInternalEnergy(rho, e) == 0.0);
      }
    }
  }
}

// Only run exception checking when we aren't offloading to the GPUs
SCENARIO("Inversion equivalence", "[SimpleMACAWEOS][Temperature]") {
  GIVEN("Initial density and specific internal energy") { 
    const Real rho_0 = 0.3;
    const Real sie_0 = 2.7;
    const Real P = eos.PressureFromDensityInternalEnergy(rho_0, sie_0);
    const Real T = eos.TemperatureFromDensityInternalEnergy(rho_0, sie_0);
    Real rho, sie;
    Real* lambda;
    eos.DensityEnergyFromPressureTemperature(P, T, lambda, rho, sie);
    sie = eos.InternalEnergyFromDensityPressure(rho, P);
    THEN("The inversion back to rho and sie from P and T should be identital.") {
      REQUIRE(isClose(rho, rho_0, REAL_TOL * rho_0));
      REQUIRE(isClose(sie, sie_0, REAL_TOL * sie_0));
    }
  }
}

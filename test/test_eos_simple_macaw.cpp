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
SCENARIO("Zero temperature on cold curve", "[SimpleMACAWEOS][Temperature]") {
  GIVEN("Parameters for a Simple MACAW EOS") {
    // Unit conversions
    //constexpr Real cm = 1.;
    //constexpr Real us = 1e-06;
    //constexpr Real Mbcc_per_g = 1e12;
    // Gruneisen parameters for copper
    constexpr Real A = 7.3;
    constexpr Real B = 3.9;
    constexpr Real Cvinf = 0.000389;
    constexpr Real v0 = 1. / 8.952;
    constexpr Real T0 = 150.;
    constexpr Real Gc = 0.5;
    // Create the EOS
    EOS host_eos = SimpleMACAW(A, B, Cvinf, v0, T0, Gc);
    //EOS eos = host_eos.GetOnDevice();
    Real rho = 0.5;
    for (int i = 0; i < 10; i++) {
      rho += rho + i; // cylce through a variety of densities
      THEN("The temperature evaluated on the cold curve should produce zero") {
        REQUIRE(eos.TemperatureFromDensityInternalEnergy(rho, eos.SieColdCurve(1.0 / rho)));
      }
    }
  }
}

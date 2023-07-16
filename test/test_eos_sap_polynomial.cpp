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
using singularity::Polynomial;

SCENARIO("Aluminium Polynomial EOS", "Check if eos returns expected values") {
  GIVEN("Parameters for a Polynomial EOS") {
    // Unit conversions
    constexpr Real cm = 1.;
    constexpr Real us = 1e-06;
    constexpr Real Mbcc_per_g = 1e12;

    constexpr Real rho0 = 2.7; // g/cc
    constexpr Real a0  = 0.0;
    constexpr Real a1  = 0.73 * cm / us;
    constexpr Real a2c = 0.73 * cm / us;
    constexpr Real a2e = 0.73 * cm / us;
    constexpr Real a3  = 0.73 * cm / us;
    constexpr Real b0  = 0.0;
    constexpr Real b1  = 0.0;
    constexpr Real b2c = 0.0;
    constexpr Real b2e = 0.0;
    constexpr Real b3  = 0.0;
    // Create the EOS
    EOS host_eos = Polynomial(rho0, a0, a1, a2c, a2e, a3, b0, b1, b2c, b2e, b3);
    EOS eos = host_eos.GetOnDevice();

    host_eos.PrintParams();
    // Calculate some stuff

    // Compare with other stuff
  }
}

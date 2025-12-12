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

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>
#include <ports-of-call/portable_errors.hpp>
#include <singularity-eos/base/indexable_types.hpp>
#include <singularity-eos/eos/eos.hpp>

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#endif

#include <test/eos_unit_test_helpers.hpp>

using singularity::DavisProducts;
using singularity::DavisReactants;
using EOS = singularity::Variant<DavisReactants, DavisProducts>;

constexpr int N = 100;
SCENARIO("Davis reactants EOS basic calls", "[DavisReactants]") {
  GIVEN("A davis reactants equation of state") {
    constexpr Real rho0 = 1.890;
    constexpr Real e0 = 4.115e10;
    constexpr Real P0 = 1.0e6;
    constexpr Real T0 = 297.0;
    constexpr Real A = 1.8e5;
    constexpr Real B = 4.6;
    constexpr Real C = 0.34;
    constexpr Real G0 = 0.56;
    constexpr Real Z = 0.0;
    constexpr Real alpha = 0.4265;
    constexpr Real Cv0 = 0.001074e10;
    EOS eos_host = DavisReactants(rho0, e0, P0, T0, A, B, C, G0, Z, alpha, Cv0);
    auto eos = eos_host.GetOnDevice();

    WHEN("The EOS is in mild expansion") {
      THEN("The Gruneisen coefficient is G0") {
        int nwrong = 0;
        portableReduce(
            "Check gruneisen coefficient in expansion", 0, N,
            PORTABLE_LAMBDA(const int i, int &nw) {
              Real rho = 0.5 * rho0;
              Real T = 10 * i;
              Real gamma = eos.GruneisenParamFromDensityTemperature(rho, T);
              nw += !(isClose(G0, gamma, 1e-12));
            },
            nwrong);
        REQUIRE(nwrong == 0);
      }
    }

    WHEN("The EOS is in compression") {
      THEN("The Gruneisen coefficient approaches G0 + 1") {
        int nwrong = 0;
        portableReduce(
            "Check gruneisen coefficient in compression", 0, N,
            PORTABLE_LAMBDA(const int i, int &nw) {
              Real rho = 10 * rho0;
              Real T = 10 * i;
              Real gamma = eos.GruneisenParamFromDensityTemperature(rho, T);
              const Real y = 1 - (rho0 / rho);
              nw += !(isClose(G0 + Z * y, gamma, 1e-12));
            },
            nwrong);
        REQUIRE(nwrong == 0);
      }
    }

    WHEN("The EOS is in extreme expansion") {
      THEN("It still produces a sensible answer") {
        int nwrong = 0;
        portableReduce(
            "Check properties in extreme expansion", 0, N,
            PORTABLE_LAMBDA(const int i, int &nw) {
              Real T = 10 * i;
              Real rho = 0;
              Real P = eos.PressureFromDensityTemperature(rho, T);
              Real E = eos.InternalEnergyFromDensityTemperature(rho, T);
              Real bmod = eos.BulkModulusFromDensityInternalEnergy(rho, T);
              nw += std::isnan(P);
              nw += std::isnan(E);
              nw += std::isnan(bmod);
              nw += (bmod < 0);
            },
            nwrong);
        REQUIRE(nwrong == 0);
      }
    }

    eos_host.Finalize();
    eos.Finalize();
  }
}

SCENARIO("Davis products EOS basic calls", "[DavisProducts]") {
  GIVEN("A davis reactants equation of state") {
    constexpr Real a = 0.798311;
    constexpr Real b = 0.58;
    constexpr Real k = 1.35;
    constexpr Real n = 2.66182;
    constexpr Real vc = 0.75419;
    constexpr Real pc = 3.2e10;
    constexpr Real Cv = 0.001072e10;
    EOS eos_host = DavisProducts(a, b, k, n, vc, pc, Cv);
    auto eos = eos_host.GetOnDevice();

    WHEN("The EOS is in extreme expansion") {
      THEN("It still produces a sensible answer") {
        int nwrong = 0;
        portableReduce(
            "Check properties in extreme expansion", 0, N,
            PORTABLE_LAMBDA(const int i, int &nw) {
              Real T = 10 * i;
              Real rho = 0;
              Real P = eos.PressureFromDensityTemperature(rho, T);
              Real E = eos.InternalEnergyFromDensityTemperature(rho, T);
              Real bmod = eos.BulkModulusFromDensityInternalEnergy(rho, T);
              nw += std::isnan(P);
              nw += std::isnan(E);
              nw += std::isnan(bmod);
              nw += (bmod < 0);
            },
            nwrong);
        REQUIRE(nwrong == 0);
      }
    }

    eos_host.Finalize();
    eos.Finalize();
  }
}

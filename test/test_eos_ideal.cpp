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

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>
#include <ports-of-call/portable_errors.hpp>
#include <singularity-eos/eos/eos.hpp>

#ifdef SINGULARITY_BUILD_CLOSURE
#include <singularity-eos/eos/singularity_eos.hpp>
#endif

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include "catch2/catch.hpp"
#endif

#include <test/eos_unit_test_helpers.hpp>

using singularity::EOS;
using singularity::IdealGas;

SCENARIO("Ideal gas entropy", "[IdealGas][Entropy]") {
  GIVEN("Parameters for an ideal gas with entropy reference states") {
    // Create ideal gas EOS ojbect
    constexpr Real Cv = 5.0;
    constexpr Real gm1 = 0.4;
    constexpr Real EntropyT0 = 100;
    constexpr Real EntropyRho0 = 1e-03;
    EOS host_eos = IdealGas(gm1, Cv, EntropyT0, EntropyRho0);
    THEN("The entropy at the reference state should be zero") {
      auto entropy = host_eos.EntropyFromDensityTemperature(EntropyRho0, EntropyT0);
      INFO("Entropy should be zero but it is " << entropy);
      CHECK(isClose(entropy, 0.0, 1.e-14));
    }
    GIVEN("A state at the reference temperature and a density whose cube root is the "
          "reference density") {
      constexpr Real T = EntropyT0;
      constexpr Real rho = 0.1; // rho**3 = EntropyRho0
      THEN("The entropy should be 2. / 3. * gm1 * Cv * log(EntropyRho0)") {
        const Real entropy_true = 2. / 3. * gm1 * Cv * log(EntropyRho0);
        auto entropy = host_eos.EntropyFromDensityTemperature(rho, T);
        INFO("Entropy: " << entropy << "  True entropy: " << entropy_true);
        CHECK(isClose(entropy, entropy_true, 1e-12));
      }
    }
    GIVEN("A state at the reference density and a temperature whose square is the "
          "reference temperature") {
      constexpr Real T = 10; // T**2 = EntropyT0
      constexpr Real rho = EntropyRho0;
      THEN("The entropy should be -1. / 2. * Cv * log(EntropyT0)") {
        const Real entropy_true = -1. / 2. * Cv * log(EntropyT0);
        auto entropy = host_eos.EntropyFromDensityTemperature(rho, T);
        INFO("Entropy: " << entropy << "  True entropy: " << entropy_true);
        CHECK(isClose(entropy, entropy_true, 1e-12));
      }
    }
  }
}

// A non-standard pattern where we do a reduction
class CheckPofRE {
 public:
  CheckPofRE(Real *P, Real *rho, Real *sie, int N) : P_(P), rho_(rho), sie_(sie), N_(N) {}
  template <typename T>
  void operator()(const T &eos) {
    Real *P = P_;
    Real *rho = rho_;
    Real *sie = sie_;
    portableReduce(
        "MyCheckPofRE", 0, N_,
        PORTABLE_LAMBDA(const int i, int &nw) {
          nw += !(isClose(P[i],
                          eos.PressureFromDensityInternalEnergy(rho[i], sie[i], nullptr),
                          1e-15));
        },
        nwrong);
  }

  // must be initialized to zero because portableReduce is a simple for loop on host
  int nwrong = 0;

 private:
  int N_;
  Real *P_;
  Real *rho_;
  Real *sie_;
};
SCENARIO("Ideal gas vector Evaluate call", "[IdealGas][Evaluate]") {
  GIVEN("An ideal gas, and some device memory") {
    constexpr Real Cv = 2.0;
    constexpr Real gm1 = 0.5;
    constexpr int N = 100;
    constexpr Real rhomin = 1;
    constexpr Real rhomax = 11;
    constexpr Real drho = (rhomax - rhomin) / (N - 1);
    constexpr Real siemin = 1;
    constexpr Real siemax = 11;
    constexpr Real dsie = (siemax - siemin) / (N - 1);

    EOS eos_host = IdealGas(gm1, Cv);
    auto eos_device = eos_host.GetOnDevice();

    Real *P = (Real *)PORTABLE_MALLOC(N * sizeof(Real));
    Real *rho = (Real *)PORTABLE_MALLOC(N * sizeof(Real));
    Real *sie = (Real *)PORTABLE_MALLOC(N * sizeof(Real));

    WHEN("We set P, rho, sie via the scalar API") {
      portableFor(
          "Set rho and sie", 0, N, PORTABLE_LAMBDA(const int i) {
            rho[i] = rhomin + drho * i; // just some diagonal in the 2d rho/sie space
            sie[i] = siemin + dsie * i;
            P[i] = eos_device.PressureFromDensityInternalEnergy(rho[i], sie[i]);
          });
      THEN("The vector Evaluate API can be used to compare") {
        CheckPofRE my_op(P, rho, sie, N);
        eos_device.Evaluate(my_op);
        REQUIRE(my_op.nwrong == 0);
      }
    }

    eos_device.Finalize();
    PORTABLE_FREE(P);
    PORTABLE_FREE(rho);
    PORTABLE_FREE(sie);
  }
}

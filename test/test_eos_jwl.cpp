//------------------------------------------------------------------------------
// © 2026. Triad National Security, LLC. All rights reserved.  This
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

// This file was made in part with generative AI.

#include <cmath>

#include <ports-of-call/portability.hpp>
#include <singularity-eos/eos/eos.hpp>

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#endif

#include <test/eos_unit_test_helpers.hpp>

using singularity::JWL;
using EOS = singularity::Variant<JWL>;

namespace {

constexpr Real A = 6.3207e14;
constexpr Real B = -4.472e10;
constexpr Real R1 = 11.3;
constexpr Real R2 = 1.13;
constexpr Real w = 0.8938;
constexpr Real rho0 = 1.895;
constexpr Real Cv = 1.31240105540897101e7;

EOS MakeJWL() { return JWL(A, B, R1, R2, w, rho0, Cv); }

} // namespace

SCENARIO("JWL EOS basic calls", "[JWL]") {
  GIVEN("A JWL EOS for detonation products") {
    EOS eos_host = MakeJWL();
    auto eos = eos_host.GetOnDevice();

    WHEN("Thermodynamic properties are evaluated at representative states") {
      constexpr int num = 4;
      int nwrong = 0;
      portableReduce(
          "Check JWL thermodynamic properties", 0, num,
          PORTABLE_LAMBDA(const int i, int &nw) {
            Real rho, temp, sie_true, press_true, bmod_true;
            if (i == 0) {
              rho = 1.42125000000000012e0;
              temp = 5.00000000000000000e2;
              sie_true = 1.94155962893205547e9;
              press_true = -1.39539537633559513e9;
              bmod_true = 3.57742905528225708e9;
            } else if (i == 1) {
              rho = 1.89500000000000002e0;
              temp = 1.00000000000000000e3;
              sie_true = 6.74299306354014778e9;
              press_true = 1.56033130109453297e10;
              bmod_true = 1.14145142327138275e11;
            } else if (i == 2) {
              rho = 2.84250000000000025e0;
              temp = 2.00000000000000000e3;
              sie_true = 3.22063102789073410e10;
              press_true = 3.83759661376449402e11;
              bmod_true = 2.65765485356934912e12;
            } else {
              rho = 3.79000000000000004e0;
              temp = 3.00000000000000000e3;
              sie_true = 1.31330196851926880e11;
              press_true = 2.33127259815571191e12;
              bmod_true = 1.27999609124014121e13;
            }

            const Real sie = eos.InternalEnergyFromDensityTemperature(rho, temp);
            nw += !isClose(sie, sie_true, 1e-12);
            nw +=
                !isClose(eos.TemperatureFromDensityInternalEnergy(rho, sie), temp, 1e-12);
            nw += !isClose(eos.PressureFromDensityTemperature(rho, temp), press_true,
                           1e-12);
            nw += !isClose(eos.PressureFromDensityInternalEnergy(rho, sie), press_true,
                           1e-12);
            nw += !isClose(eos.SpecificHeatFromDensityTemperature(rho, temp), Cv, 1e-12);
            nw +=
                !isClose(eos.SpecificHeatFromDensityInternalEnergy(rho, sie), Cv, 1e-12);
            nw += !isClose(eos.BulkModulusFromDensityTemperature(rho, temp), bmod_true,
                           1e-12);
            nw += !isClose(eos.BulkModulusFromDensityInternalEnergy(rho, sie), bmod_true,
                           1e-12);
            nw += !isClose(eos.GruneisenParamFromDensityTemperature(rho, temp), w, 1e-12);
            nw +=
                !isClose(eos.GruneisenParamFromDensityInternalEnergy(rho, sie), w, 1e-12);
          },
          nwrong);

      THEN("The values agree with the reference results") { REQUIRE(nwrong == 0); }
    }

    WHEN("Density and energy are recovered from pressure and temperature") {
      constexpr int num = 16;
      int nwrong = 0;
      portableReduce(
          "Check JWL pressure-temperature inversion", 0, num,
          PORTABLE_LAMBDA(const int i, int &nw) {
            const Real rho = rho0 * (1.0 + static_cast<Real>(i) / num);
            const Real temp = 500.0 + 100.0 * i;
            nw += !CheckRhoSieFromPT(eos, rho, temp);
          },
          nwrong);

      THEN("The original density and energy are recovered") { REQUIRE(nwrong == 0); }
    }

    eos.Finalize();
    eos_host.Finalize();
  }
}

SCENARIO("JWL isentropic bulk modulus agrees with finite differences",
         "[JWL][BulkModulus]") {
  GIVEN("A JWL EOS and a range of physical states") {
    EOS eos_host = MakeJWL();
    auto eos = eos_host.GetOnDevice();

    constexpr int num = 32;
    int nwrong = 0;
    portableReduce(
        "Compare JWL bulk modulus to finite differences", 0, num,
        PORTABLE_LAMBDA(const int i, int &nw) {
          const Real rho = rho0 * (0.8 + 0.04 * i);
          const Real temp = 500.0 + 75.0 * i;
          const Real sie = eos.InternalEnergyFromDensityTemperature(rho, temp);
          const Real press = eos.PressureFromDensityInternalEnergy(rho, sie);
          const Real drho = 1e-6 * rho;
          const Real dsie = 1e-6 * std::abs(sie);

          const Real dpdr_e = (eos.PressureFromDensityInternalEnergy(rho + drho, sie) -
                               eos.PressureFromDensityInternalEnergy(rho - drho, sie)) /
                              (2.0 * drho);
          const Real dpde_r = (eos.PressureFromDensityInternalEnergy(rho, sie + dsie) -
                               eos.PressureFromDensityInternalEnergy(rho, sie - dsie)) /
                              (2.0 * dsie);
          const Real bmod_fd = rho * dpdr_e + press / rho * dpde_r;
          const Real bmod = eos.BulkModulusFromDensityInternalEnergy(rho, sie);
          nw += !isClose(bmod, bmod_fd, 1e-8);
        },
        nwrong);

    THEN("The analytic and finite-difference results agree") { REQUIRE(nwrong == 0); }

    eos.Finalize();
    eos_host.Finalize();
  }
}

SCENARIO("JWL internal energy from density and pressure", "[JWL][SieFromRhoP]") {
  GIVEN("A state from Lee et al. (2013). 10.1016/j.jcp.2013.03.046") {
    // The source uses g, cm, and microseconds, so pressures in Mbar and specific
    // energies in (cm/us)^2 both acquire a factor of 1e12 when converted to cgs.
    constexpr Real paper_units_to_cgs = 1e12;
    constexpr Real A = 632.1 * paper_units_to_cgs;
    constexpr Real B = -0.04472 * paper_units_to_cgs;
    constexpr Real R1 = 11.3;
    constexpr Real R2 = 1.13;
    constexpr Real w = 0.8938;
    constexpr Real rho0 = 1.905;
    constexpr Real Cv = 0.00003888 * paper_units_to_cgs;
    constexpr Real rho = 3.810;
    constexpr Real press = 2.0 * paper_units_to_cgs;
    constexpr Real sie_true = 0.0333355 * paper_units_to_cgs;

    EOS eos_host = JWL(A, B, R1, R2, w, rho0, Cv);
    auto eos = eos_host.GetOnDevice();

    WHEN("Internal energy is computed from density and pressure") {
      Real sie = 0.0;
      portableReduce(
          "Compute JWL internal energy from density and pressure", 0, 1,
          PORTABLE_LAMBDA(const int, Real &result) {
            eos.InternalEnergyFromDensityPressure(rho, press, result);
          },
          sie);

      THEN("It agrees with the published state") {
        INFO("Computed sie: " << sie << ", expected sie: " << sie_true);
        REQUIRE(isClose(sie, sie_true, 1e-6));
      }
    }

    eos.Finalize();
    eos_host.Finalize();
  }
}

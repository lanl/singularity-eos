//------------------------------------------------------------------------------
// © 2024. Triad National Security, LLC. All rights reserved.  This
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

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#endif

#include <test/eos_unit_test_helpers.hpp>

using singularity::IdealElectrons;
using singularity::IdealGas;
using singularity::MeanAtomicProperties;
using singularity::ZSplitE;
using singularity::ZSplitI;

using EOS =
    singularity::Variant<IdealGas, IdealElectrons, ZSplitI<IdealGas>, ZSplitE<IdealGas>>;

SCENARIO("ZSplit of Ideal Gas", "[ZSplit][IdealGas][IdealElectrons]") {
  GIVEN("An ideal gas EOS") {
    constexpr Real gm1 = (5. / 3.) - 1.;

    constexpr Real Abar = 26; // Iron
    constexpr Real Zbar = 55.845;
    const MeanAtomicProperties azbar(Abar, Zbar);

    constexpr Real Cv =
        (Zbar + 1) * IdealElectrons::kb / (gm1 * Abar * IdealElectrons::mp);

    EOS eos_t_h = IdealGas(gm1, Cv, azbar);
    EOS eos_ze_h = eos_t_h.Modify<ZSplitE>();
    EOS eos_zi_h = eos_t_h.Modify<ZSplitI>();
    EOS eos_ie_h = IdealElectrons(azbar);

    auto eos_t = eos_t_h.GetOnDevice();
    auto eos_ze = eos_ze_h.GetOnDevice();
    auto eos_zi = eos_zi_h.GetOnDevice();
    auto eos_ie = eos_ie_h.GetOnDevice();
    WHEN("We evaluate at fixed density and temperature for several ionization states") {
      constexpr Real rho = 1.0;
      constexpr Real temp = 1e3;
      constexpr Real zmin = 0;
      constexpr Real zmax = Zbar;
      constexpr std::size_t NZ = 1000;
      constexpr Real dz = (zmax - zmin) / (NZ - 1);

      int nwrong = 0;
      portableReduce(
          "Zsplit check P and sie", 0, NZ,
          PORTABLE_LAMBDA(const int i, int &nw) {
            Real Z = zmin + dz * i;
            Real lambda[1] = {Z};

            Real Pt = eos_t.PressureFromDensityTemperature(rho, temp, lambda);
            Real Pi = eos_zi.PressureFromDensityTemperature(rho, temp, lambda);
            Real Pe = eos_ze.PressureFromDensityTemperature(rho, temp, lambda);
            if (!isClose(Pt, Pi + Pe, 1e-12)) {
              printf("Bad P Sum! %d %.14e %.14e %.14e\n", i, Pt, Pi, Pe);
              nw += 1;
            }

            Real et = eos_t.InternalEnergyFromDensityTemperature(rho, temp, lambda);
            Real ei = eos_zi.InternalEnergyFromDensityTemperature(rho, temp, lambda);
            Real ee = eos_ze.InternalEnergyFromDensityTemperature(rho, temp, lambda);
            if (!isClose(et, ei + ee, 1e-12)) {
              printf("Bad sie Sum! %d %.14e %.14e %.14e\n", i, et, ei, ee);
              nw += 1;
            }
          },
          nwrong);
      THEN("Pe + Pi = Pt") { REQUIRE(nwrong == 0); }

      WHEN("Ideal electrons are fully ionized") {
        constexpr Real rho = 1.0;
        constexpr Real Tmin = 3e3;
        constexpr Real Tmax = 3e5;
        constexpr int NT = 1000;
        constexpr Real dt = (Tmax - Tmin) / (NT - 1);

        int nwrong = 0;
        portableReduce(
            "Ideal electrons check full ionization", 0, NT,
            PORTABLE_LAMBDA(const int i, int &nw) {
              const Real T = Tmin + i * dt;
              Real lambda[1] = {Zbar};

              Real P_fi = eos_ze.PressureFromDensityTemperature(rho, T, lambda);
              Real P_ie = eos_ie.PressureFromDensityTemperature(rho, T, lambda);
              if (!isClose(P_fi, P_ie, 1e-12)) {
                printf("Bad P ionized! %d %.14e %.14e\n", i, P_fi, P_ie);
                nw += 1;
              }

              Real e_fi = eos_ze.InternalEnergyFromDensityTemperature(rho, T, lambda);
              Real e_ie = eos_ie.InternalEnergyFromDensityTemperature(rho, T, lambda);
              if (!isClose(e_fi, e_ie, 1e-12)) {
                printf("Bad P ionized! %d %.14e %.14e\n", i, e_fi, e_ie);
                nw += 1;
              }
            },
            nwrong);
        THEN("Pe_fi == PE_ie") { REQUIRE(nwrong == 0); }
      }

      eos_t.Finalize();
      eos_ze.Finalize();
      eos_zi.Finalize();
      eos_ie.Finalize();
    }
  }
}

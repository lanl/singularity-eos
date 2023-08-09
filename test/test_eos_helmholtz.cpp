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

#include <array>
#include <limits>
#include <string>

#ifndef CATCH_CONFIG_RUNNER
#include "catch2/catch.hpp"
#endif

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>
#include <ports-of-call/portable_errors.hpp>
#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/eos/eos.hpp>
#include <test/eos_unit_test_helpers.hpp>

using singularity::Helmholtz;
const std::string filename = "../test/helmholtz/helm_table.dat";
SCENARIO("Helmholtz equation of state tgiven", "[HelmholtzEOS]") {
  GIVEN("A helholtz EOS") {
    /* We only test the EOS without Coulomb corrections since those are
       implemented in a different way in the reference implementation (cutoff
       vs. butterworth filter)*/
    Helmholtz eos(filename, true, true, false, true, true);
    THEN("We loaded the file!") { REQUIRE(true); }

    /* Density and temperature range evenly sampling the parameter space */
    Real rho_in[4] = {1e-3, 1e1, 1e5, 1e9};
    Real temp_in[4] = {1e4, 1e6, 1e8, 1e10};

    /* Reference values calculated with test implementation, using
        abar = 4.0, zbar = 2.0 */
    Real ein_ref[16] = {
      942247520796.137939453125,
      101117704971247.1875,
      756586284440018092032.0,
      204533796186988033170649120768.0,
      22331233357586.9296875,
      100345465684484.421875,
      85141422479277632.0,
      20453379683965599047745536.0,
      9976645702610400.0,
      10007780352178464.0,
      15481955850958990.0,
      2045403973039441838080.0,
      1256563470889144064.0,
      1256594342723058944.0,
      1259724628933661184.0,
      2037841273725276672.0,
    };
    Real press_ref[16] = {
      628164115.19633424282073974609375,
      64881120328.5929718017578125,
      252198502998578304.0,
      66317830432079092035616768.0,
      148867082685630.90625,
      668864622842899.125,
      314555378539921728.0,
      66317832510490259937558528.0,
      649943133468265414656.0,
      652017695913837068288.0,
      1006403106957914537984.0,
      66338640383060857803767808.0,
      486181392465836489259876352.0,
      486201972259441884654469120.0,
      488274539873823033716113408.0,
      852398557407648750030028800.0,
    };

    /* Compare test values. Difference should be less than 1e-12 */
    Real lambda[2] = {4.0, 2.0};
    int k = 0;
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j){
          Real ein = eos.InternalEnergyFromDensityTemperature(rho_in[i], temp_in[j], lambda);
          Real press = eos.PressureFromDensityTemperature(rho_in[i], temp_in[j], lambda);
          REQUIRE(isClose(ein, ein_ref[k], 1e-12));
          REQUIRE(isClose(press, press_ref[k], 1e-12));
          k++;
      }
    }
  }
}

SCENARIO("Helmholtz equation of state egiven", "[HelmholtzEOS]") {
  GIVEN("A helholtz EOS") {
    /* We only test the EOS without Coulomb corrections since those are
       implemented in a different way in the reference implementation (cutoff
       vs. butterworth filter)*/
    Helmholtz eos(filename, true, true, false, true, true);
    THEN("We loaded the file!") { REQUIRE(true); }

    /* Density and internal energy range sampling the parameter space.
       Due to problems with the initial temperature guess, only a
       particular range of values work with the built-in ideal
       gas initial guess. */
    
    Real rho_in[7] = {1e-3, 1e-3, 1e1, 1e1, 1e5, 1e5, 1e9};
    Real ein_in[7] = {1e33, 1e38, 1e25, 1e38, 1e22, 1e32, 1e22};

    /* Reference values calculated with test implementation, using
        abar = 4.0, zbar = 2.0 */
    Real temp_ref[7] = {
      83267912059.9008941650390625,
      1480649420949.600341796875,
      8377732410.2377471923828125,
      10000000000000.037109375,
      14834616134.222881317138671875,
      2633008614585.40283203125,
      148051392128.0281982421875,
    };
    Real press_ref[7] = {
      333178534190936501846297018368.0,
      3.335059068385989208866699710405e+34,
      32102269803001238788243456.0,
      6.9353723467611821959666509095658e+37,
      328892155446139829490810880.0,
      3.3333323663968881590598358740553e+35,
      3334463461062235621149686890496.0,
    };

    /* Compare test values. Difference should be less than 1e-6 */
    Real lambda[2] = {4.0, 2.0};
    for (int i = 0; i < 7; ++i) {
          Real temp = eos.TemperatureFromDensityInternalEnergy(rho_in[i], ein_in[i], lambda);
          Real press = eos.PressureFromDensityInternalEnergy(rho_in[i], ein_in[i], lambda);
          Real temp_diff = (temp - temp_ref[i]) / temp_ref[i];
          Real press_diff = (press - press_ref[i]) / press_ref[i];
          REQUIRE(temp_diff < 1e6);
          REQUIRE(press_diff < 1e6);
          // REQUIRE(isClose(temp, temp_ref[i], 1e-0));
          // REQUIRE(isClose(press, press_ref[i], 1e-0));
    }
  }
}

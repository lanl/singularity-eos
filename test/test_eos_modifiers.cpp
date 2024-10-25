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

#include <limits>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>
#include <ports-of-call/portable_errors.hpp>
#include <singularity-eos/base/fast-math/logs.hpp>
#include <singularity-eos/base/root-finding-1d/root_finding.hpp>
#include <singularity-eos/base/serialization_utils.hpp>
#include <singularity-eos/eos/eos.hpp>
#include <singularity-eos/eos/eos_builder.hpp>

#ifdef SINGULARITY_BUILD_CLOSURE
#include <singularity-eos/eos/singularity_eos.hpp>
#endif

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#endif

#include <test/eos_unit_test_helpers.hpp>

namespace eos_units_init = singularity::eos_units_init;

namespace EOSBuilder = singularity::EOSBuilder;
namespace thermalqs = singularity::thermalqs;
namespace variadic_utils = singularity::variadic_utils;

using EOSBuilder::Modify;
using singularity::BilinearRampEOS;
using singularity::IdealGas;
using singularity::RelativisticEOS;
using singularity::ScaledEOS;
using singularity::ShiftedEOS;
using singularity::UnitSystem;

#ifndef SINGULARITY_BUILD_CLOSURE
// recreate variadic list
template <typename... Ts>
using tl = variadic_utils::type_list<Ts...>;

template <template <typename> class... Ts>
using al = variadic_utils::adapt_list<Ts...>;

// transform variadic list: applies modifiers to eos's
using variadic_utils::transform_variadic_list;

static constexpr const auto full_eos_list = tl<IdealGas>{};
static constexpr const auto relativistic_eos_list = tl<IdealGas>{};
static constexpr const auto unit_system_eos_list = tl<IdealGas>{};
static constexpr const auto apply_to_all = al<ScaledEOS, ShiftedEOS>{};
static constexpr const auto unit_system =
    transform_variadic_list(unit_system_eos_list, al<UnitSystem>{});
// variadic list of eos's with shifted or scaled modifiers
static constexpr const auto shifted_1 =
    transform_variadic_list(full_eos_list, al<ShiftedEOS>{});
static constexpr const auto scaled_1 =
    transform_variadic_list(full_eos_list, al<ScaledEOS>{});
// variadic list of Relativistic<T>'s
static constexpr const auto relativistic =
    transform_variadic_list(relativistic_eos_list, al<RelativisticEOS>{});
// relativistic and unit system modifiers
static constexpr const auto unit_or_rel =
    variadic_utils::concat(unit_system, relativistic);
// variadic list of eos with shifted, relativistic or unit system modifiers
static constexpr const auto shifted_of_unit_or_rel =
    transform_variadic_list(unit_or_rel, al<ShiftedEOS>{});
// combined list of all shifted EOS
static constexpr const auto shifted =
    variadic_utils::concat(shifted_1, shifted_of_unit_or_rel);
// variadic list of eos with scaled, relativistic or unit system modifiers
static constexpr const auto scaled_of_unit_or_rel =
    transform_variadic_list(unit_or_rel, al<ScaledEOS>{});
// variadic list of Scaled<Shifted<T>>'s
static constexpr const auto scaled_of_shifted =
    transform_variadic_list(shifted, al<ScaledEOS>{});
// combined list of all scaled EOS
static constexpr const auto scaled =
    variadic_utils::concat(scaled_1, scaled_of_unit_or_rel, scaled_of_shifted);
// create combined list
static constexpr const auto combined_list_1 =
    variadic_utils::concat(full_eos_list, shifted, scaled, unit_or_rel);
// make a ramped eos of everything
static constexpr const auto ramped_all =
    transform_variadic_list(combined_list_1, al<BilinearRampEOS>{});
// final combined list
static constexpr const auto combined_list =
    variadic_utils::concat(combined_list_1, ramped_all);
using EOS = typename decltype(tl_to_Variant(combined_list))::vt;
#endif

SCENARIO("EOS Builder and Modifiers", "[EOSBuilder][Modifiers][IdealGas]") {

  GIVEN("Parameters for a shifted and scaled ideal gas") {
    constexpr Real Cv = 2.0;
    constexpr Real gm1 = 0.5;
    constexpr Real scale = 2.0;
    constexpr Real shift = 0.1;
    constexpr Real rho = 2.0;
    constexpr Real sie = 0.5;
    WHEN("We use the EOSBuilder") {
      EOS eos = IdealGas(gm1, Cv);
      eos = eos.Modify<ShiftedEOS>(shift);
      eos = eos.Modify<ScaledEOS>(scale);
      THEN("The shift and scale parameters pass through correctly") {
        REQUIRE(eos.PressureFromDensityInternalEnergy(rho, sie) == 0.3);
      }
      THEN("We can UnmodifyOnce to get the shifted EOS object") {
        EOS shifted = eos.UnmodifyOnce();
        REQUIRE(shifted.IsType<ShiftedEOS<IdealGas>>());
        AND_THEN("We can extract the unmodified object") {
          EOS unmod = eos.GetUnmodifiedObject();
          REQUIRE(unmod.IsType<IdealGas>());
        }
      }
      WHEN("We add a ramp") {
        // EOSBuilder::params_t ramp_params;
        Real r0 = 1;
        Real a = 1;
        Real b = 0;
        Real c = 0;
        THEN("The EOS is constructed correctly") {
          auto eos_ramped = Modify<BilinearRampEOS>(eos, r0, a, b, c);
        }
      }
    }
#ifdef SINGULARITY_BUILD_CLOSURE
    WHEN("We construct a non-modifying modifier") {
      EOS ig = IdealGas(gm1, Cv);
      EOS igsh = ScaledEOS<IdealGas>(IdealGas(gm1, Cv), 1.0);
      EOS igsc = ShiftedEOS<IdealGas>(IdealGas(gm1, Cv), 0.0);
      EOS igra;
      // test out the c interface
      int enabled[4] = {0, 0, 1, 0};
      Real vals[6] = {0.0, 0.0, 1.e9, 1.0, 2.0, 1.0};
      init_sg_IdealGas(0, &igra, gm1, Cv, enabled, vals);
      THEN("The modified EOS should produce equivalent results") {
        compare_two_eoss(igsh, ig);
        compare_two_eoss(igsc, ig);
        compare_two_eoss(igra, ig);
      }
    }
    WHEN("We construct a ramp from a p-alpha model") {
      const Real Pe = 5.e7, Pc = 1.e8;
      const Real alpha0 = 1.5;
      const Real T0 = 293.0;
      int enabled[4] = {0, 0, 0, 1};
      Real vals[6] = {0.0, 0.0, alpha0, Pe, Pc, 0.0};
      const Real rho0 = 1.e6 / (gm1 * Cv * T0);
      EOS igra;
      const Real r0 = rho0 / alpha0;
      const Real r1 = Pc / (gm1 * Cv * T0);
      const Real rmid = Pe / (gm1 * Cv * T0 * alpha0);
      // P(alpha0 * rmid)
      const Real P_armid = alpha0 * gm1 * Cv * rmid * T0;
      init_sg_IdealGas(0, &igra, gm1, Cv, enabled, vals);
      // construct ramp params and evaluate directly for test
      const Real a = r0 * Pe / (rmid - r0);
      const Real b = r0 * (Pc - Pe) / (r1 - rmid);
      const Real c = (Pc * rmid - Pe * r1) / (r0 * (Pc - Pe));
      // density in the middle of the first slope
      const Real rho_t1 = 0.5 * (r0 + rmid);
      // density in the middle of the second slope
      const Real rho_t2 = 0.5 * (rmid + r1);
      // P (rho_t1) note that r0 = rho0 / alpha0
      const Real Prhot1 = a * (rho_t1 / r0 - 1.0);
      // P (rho_t2)
      const Real Prhot2 = b * (rho_t2 / r0 - c);
      // bmod (rho_t1)
      const Real bmodrt1 = rho_t1 * a / r0;
      // bmod (rho_t2)
      const Real bmodrt2 = rho_t2 * b / r0;
      THEN("P_eos(alpha_0*rmid, T0) = P_ramp(rmid,T0)") {
        INFO("P_eos(alpha_0*rmid, T0): "
             << P_armid
             << " P_ramp(rmid, T0): " << igra.PressureFromDensityTemperature(rmid, T0));
        REQUIRE(isClose(P_armid, igra.PressureFromDensityTemperature(rmid, T0), 1.e-12));
      }
      THEN("We obtain correct ramp behavior in P(rho) for rho <r0, [r0,rmid], [rmid,r1] "
           "and >r1") {
        // also check pressures on ramp
        INFO("reference P((r0+rmid)/2, T0): "
             << Prhot1 << " test P((r0+rmid)/2, T0): "
             << igra.PressureFromDensityTemperature(rho_t1, T0));
        REQUIRE(isClose(Prhot1, igra.PressureFromDensityTemperature(rho_t1, T0), 1.e-12));
        INFO("reference P((rmid+r1)/2, T0): "
             << Prhot2 << " test P((rmid+r1)/2, T0): "
             << igra.PressureFromDensityTemperature(rho_t2, T0));
        REQUIRE(isClose(Prhot2, igra.PressureFromDensityTemperature(rho_t2, T0), 1.e-12));
        // check pressure below and beyond ramp matches unmodified ideal gas
        INFO("reference P(0.8*r0, T0): "
             << 0.8 * r0 * gm1 * Cv * T0 << " test P(0.8*r0, T0): "
             << igra.PressureFromDensityTemperature(0.8 * r0, T0));
        REQUIRE(isClose(0.8 * r0 * gm1 * Cv * T0,
                        igra.PressureFromDensityTemperature(0.8 * r0, T0), 1.e-12));
        INFO("reference P(1.2*r1, T0): "
             << 1.2 * r1 * gm1 * Cv * T0 << " test P(1.2*r1, T0): "
             << igra.PressureFromDensityTemperature(1.2 * r1, T0));
        REQUIRE(isClose(1.2 * r1 * gm1 * Cv * T0,
                        igra.PressureFromDensityTemperature(1.2 * r1, T0), 1.e-12));
      }
      THEN("We obtain correct ramp behavior in bmod(rho) for rho <r0, [r0,rmid], "
           "[rmid,r1] and >r1") {
        // check bulk moduli on both pieces of ramp
        INFO("reference bmod((r0+rmid)/2, T0): "
             << bmodrt1 << " test bmod((r0+rmid)/2, T0): "
             << igra.BulkModulusFromDensityTemperature(rho_t1, T0));
        REQUIRE(
            isClose(bmodrt1, igra.BulkModulusFromDensityTemperature(rho_t1, T0), 1.e-12));
        INFO("reference bmod((rmid+r1)/2, T0): "
             << bmodrt2 << " test bmod((rmid+r1)/2, T0): "
             << igra.BulkModulusFromDensityTemperature(rho_t2, T0));
        REQUIRE(
            isClose(bmodrt2, igra.BulkModulusFromDensityTemperature(rho_t2, T0), 1.e-12));
        // check bulk modulus below and beyond ramp matches unmodified ideal gas
        INFO("reference bmod(0.8*r0, T0): "
             << 0.8 * r0 * gm1 * (gm1 + 1.0) * Cv * T0 << " test bmod(0.8*r0, T0): "
             << igra.BulkModulusFromDensityTemperature(0.8 * r0, T0));
        REQUIRE(isClose(0.8 * r0 * gm1 * (gm1 + 1.0) * Cv * T0,
                        igra.BulkModulusFromDensityTemperature(0.8 * r0, T0), 1.e-12));
        INFO("reference bmod(1.2*r1, T0): "
             << 1.2 * r1 * gm1 * (gm1 + 1.0) * Cv * T0 << " test bmod(1.2*r1, T0): "
             << igra.BulkModulusFromDensityTemperature(1.2 * r1, T0));
        REQUIRE(isClose(1.2 * r1 * gm1 * (gm1 + 1.0) * Cv * T0,
                        igra.BulkModulusFromDensityTemperature(1.2 * r1, T0), 1.e-12));
      }
    }
#endif // SINGULARITY_BUILD_CLOSURE
  }
}

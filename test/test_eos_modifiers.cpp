//------------------------------------------------------------------------------
// © 2021-2026. Triad National Security, LLC. All rights reserved.  This
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

// This file was generated in part with the assistance of generative AI

#include <limits>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>
#include <ports-of-call/portable_errors.hpp>
#include <singularity-eos/base/fast-math/logs.hpp>
#include <singularity-eos/base/robust_utils.hpp>
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

/* A toy version of ideal gas where bounds have been placed on it so
   we can check these are passed through modifiers properly.
 */
class BoundedGas : public IdealGas {
 public:
  BoundedGas() = default;
  PORTABLE_INLINE_FUNCTION
  BoundedGas(Real gm1, Real Cv) : IdealGas(gm1, Cv) {}

  PORTABLE_FORCEINLINE_FUNCTION
  Real MinimumDensity() const { return 1e-2; }
  PORTABLE_FORCEINLINE_FUNCTION
  Real MinimumTemperature() const { return 1e-3; }

  PORTABLE_FORCEINLINE_FUNCTION
  Real MaximumDensity() const { return 1e10; }

  PORTABLE_FORCEINLINE_FUNCTION
  Real MinimumPressure() const {
    return PressureFromDensityTemperature(MinimumDensity(), MinimumTemperature());
  }
  // Gruneisen EOS's often have a maximum density, which implies a maximum pressure.
  PORTABLE_FORCEINLINE_FUNCTION
  Real MaximumPressureAtTemperature([[maybe_unused]] const Real T) const {
    return MaximumDensity();
  }

  PORTABLE_INLINE_FUNCTION
  Real RhoPmin(const Real /*temp*/) const { return MinimumDensity(); }
};

#ifndef SINGULARITY_BUILD_CLOSURE
// recreate variadic list
template <typename... Ts>
using tl = variadic_utils::type_list<Ts...>;

template <template <typename> class... Ts>
using al = variadic_utils::adapt_list<Ts...>;

// transform variadic list: applies modifiers to eos's
using variadic_utils::transform_variadic_list;

static constexpr const auto full_eos_list = tl<IdealGas>{};
// modifiers that get applied to all eos's
static constexpr const auto apply_to_all = al<ScaledEOS, ShiftedEOS>{};
// variadic list of eos's with shifted or scaled modifiers
static constexpr const auto shifted =
    transform_variadic_list(full_eos_list, al<ShiftedEOS>{});
static constexpr const auto scaled_1 =
    transform_variadic_list(full_eos_list, al<ScaledEOS>{});
// variadic list of Scaled<Shifted<T>>'s
static constexpr const auto scaled_of_shifted =
    transform_variadic_list(shifted, al<ScaledEOS>{});
// combined list of all scaled EOS
static constexpr const auto scaled =
    singularity::variadic_utils::concat(scaled_1, scaled_of_shifted);
// create combined list
static constexpr const auto combined_list_1 =
    singularity::variadic_utils::concat(full_eos_list, shifted, scaled);
// make a ramped eos of everything
static constexpr const auto ramped_all =
    transform_variadic_list(combined_list_1, al<BilinearRampEOS>{});
// final combined list
static constexpr const auto combined_list =
    singularity::variadic_utils::concat(combined_list_1, ramped_all);

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
          eos_ramped.PrintParams();
        }
      }

      WHEN("We compute PT derivatives from preferred for the modified EOS") {
        // Choose a representative state
        constexpr Real test_rho = 1.5;
        constexpr Real test_T = 2.0;

        // Compute primitive quantities
        const Real P = eos.PressureFromDensityTemperature(test_rho, test_T);
        const Real sie_val = eos.InternalEnergyFromDensityTemperature(test_rho, test_T);

        // Derivatives from the EOS preferred interface
        Real dedP_T, drdP_T, dedT_P, drdT_P;
        eos.PTDerivativesFromPreferred(test_rho, sie_val, P, test_T,
                                       static_cast<Real *>(nullptr), dedP_T, drdP_T,
                                       dedT_P, drdT_P);

        // Finite‑difference reference
        constexpr Real derivative_eps = 3e-6;
        Real dedP_T_fd, drdP_T_fd, dedT_P_fd, drdT_P_fd;
        singularity::eos_base::PTDerivativesByFiniteDifferences(
            eos, test_rho, sie_val, P, test_T, derivative_eps,
            static_cast<Real *>(nullptr), dedP_T_fd, drdP_T_fd, dedT_P_fd, drdT_P_fd);

        // Compare each component
        THEN("dedP_T matches finite‑difference reference") {
          bool close = isClose(dedP_T, dedP_T_fd);
          if (!close) {
            printf("dedP_T mismatch: eos = %.14e, fd = %.14e, diff = %.14e\n", dedP_T,
                   dedP_T_fd, dedP_T - dedP_T_fd);
          }
          REQUIRE(close);
        }
        THEN("drdP_T matches finite‑difference reference") {
          bool close = isClose(drdP_T, drdP_T_fd);
          if (!close) {
            printf("drdP_T mismatch: eos = %.14e, fd = %.14e, diff = %.14e\n", drdP_T,
                   drdP_T_fd, drdP_T - drdP_T_fd);
          }
          REQUIRE(close);
        }
        THEN("dedT_P matches finite‑difference reference") {
          bool close = isClose(dedT_P, dedT_P_fd);
          if (!close) {
            printf("dedT_P mismatch: eos = %.14e, fd = %.14e, diff = %.14e\n", dedT_P,
                   dedT_P_fd, dedT_P - dedT_P_fd);
          }
          REQUIRE(close);
        }
        THEN("drdT_P matches finite‑difference reference") {
          bool close = isClose(drdT_P, drdT_P_fd);
          if (!close) {
            printf("drdT_P mismatch: eos = %.14e, fd = %.14e, diff = %.14e\n", drdT_P,
                   drdT_P_fd, drdT_P - drdT_P_fd);
          }
          REQUIRE(close);
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
        printf("P_eos(alpha_0*rmid, T0): %.14e  P_ramp(rmid, T0): %.14e\n", P_armid,
               igra.PressureFromDensityTemperature(rmid, T0));
        REQUIRE(isClose(P_armid, igra.PressureFromDensityTemperature(rmid, T0), 1.e-12));
      }
      THEN("We obtain correct ramp behavior in P(rho) for rho <r0, [r0,rmid], [rmid,r1] "
           "and >r1") {
        // also check pressures on ramp
        printf("reference P((r0+rmid)/2, T0): %.14e  test P((r0+rmid)/2, T0): %.14e\n",
               Prhot1, igra.PressureFromDensityTemperature(rho_t1, T0));
        REQUIRE(isClose(Prhot1, igra.PressureFromDensityTemperature(rho_t1, T0), 1.e-12));
        printf("reference P((rmid+r1)/2, T0): %.14e  test P((rmid+r1)/2, T0): %.14e\n",
               Prhot2, igra.PressureFromDensityTemperature(rho_t2, T0));
        REQUIRE(isClose(Prhot2, igra.PressureFromDensityTemperature(rho_t2, T0), 1.e-12));
        // check pressure below and beyond ramp matches unmodified ideal gas
        printf("reference P(0.8*r0, T0): %.14e  test P(0.8*r0, T0): %.14e\n",
               0.8 * r0 * gm1 * Cv * T0,
               igra.PressureFromDensityTemperature(0.8 * r0, T0));
        REQUIRE(isClose(0.8 * r0 * gm1 * Cv * T0,
                        igra.PressureFromDensityTemperature(0.8 * r0, T0), 1.e-12));
        printf("reference P(1.2*r1, T0): %.14e  test P(1.2*r1, T0): %.14e\n",
               1.2 * r0 * gm1 * Cv * T0,
               igra.PressureFromDensityTemperature(1.2 * r0, T0));
        REQUIRE(isClose(1.2 * r0 * gm1 * Cv * T0,
                        igra.PressureFromDensityTemperature(1.2 * r0, T0), 1.e-12));
      }
      THEN("We obtain correct ramp behavior in bmod(rho) for rho <r0, [r0,rmid], "
           "[rmid,r1] and >r0") {
        // check bulk moduli on both pieces of ramp
        printf(
            "reference bmod((r0+rmid)/2, T0): %.14e  test bmod((r0+rmid)/2, T0): %.14e\n",
            bmodrt1, igra.BulkModulusFromDensityTemperature(rho_t1, T0));
        REQUIRE(
            isClose(bmodrt1, igra.BulkModulusFromDensityTemperature(rho_t1, T0), 1.e-12));
        printf(
            "reference bmod((rmid+r1)/2, T0): %.14e  test bmod((rmid+r1)/2, T0): %.14e\n",
            bmodrt2, igra.BulkModulusFromDensityTemperature(rho_t2, T0));
        REQUIRE(
            isClose(bmodrt2, igra.BulkModulusFromDensityTemperature(rho_t2, T0), 1.e-12));
        // check bulk modulus below and beyond ramp matches unmodified ideal gas
        printf("reference bmod(0.8*r0, T0): %.14e  test bmod(0.8*r0, T0): %.14e\n",
               0.8 * r0 * gm1 * (gm1 + 1.0) * Cv * T0,
               igra.BulkModulusFromDensityTemperature(0.8 * r0, T0));
        REQUIRE(isClose(0.8 * r0 * gm1 * (gm1 + 1.0) * Cv * T0,
                        igra.BulkModulusFromDensityTemperature(0.8 * r0, T0), 1.e-12));
        printf("reference bmod(1.2*r1, T0): %.14e  test bmod(1.2*r1, T0): %.14e\n",
               1.2 * r1 * gm1 * (gm1 + 1.0) * Cv * T0,
               igra.BulkModulusFromDensityTemperature(1.2 * r1, T0));
        REQUIRE(isClose(1.2 * r1 * gm1 * (gm1 + 1.0) * Cv * T0,
                        igra.BulkModulusFromDensityTemperature(1.2 * r1, T0), 1.e-12));
      }
    }
#endif // SINGULARITY_BUILD_CLOSURE
  }
}

SCENARIO("Modifiers propagate introspection bounds correctly", "[Modifiers]") {
  GIVEN("A BoundedGas EOS with known bounds") {
    constexpr Real Cv = 2.0;
    constexpr Real gm1 = 0.5;
    constexpr Real shift = 0.1;
    constexpr Real scale = 2.0;
    constexpr Real EPS = 10 * singularity::robust::EPS();

    BoundedGas bg(gm1, Cv);
    const Real base_min_rho = bg.MinimumDensity();
    const Real base_max_rho = bg.MaximumDensity();
    const Real base_min_temp = bg.MinimumTemperature();
    const Real base_min_pres = bg.MinimumPressure();
    const Real base_max_pres = bg.MaximumPressureAtTemperature(0.0);
    const Real base_rho_pmin = bg.RhoPmin(0.0);

    AND_GIVEN("A shifted, scaled EOS") {
      auto eos = ScaledEOS<ShiftedEOS<BoundedGas>>(
          ShiftedEOS<BoundedGas>(BoundedGas(gm1, Cv), shift), scale);

      THEN("MinimumDensity is scaled by the scale factor") {
        REQUIRE(isClose(eos.MinimumDensity(), base_min_rho / scale, 1.e-12));
      }

      THEN("MaximumDensity is scaled by the scale factor") {
        REQUIRE(isClose(eos.MaximumDensity(), base_max_rho / scale, 1.e-12));
      }

      THEN("MinimumTemperature is unchanged by shift/scale") {
        REQUIRE(isClose(eos.MinimumTemperature(), base_min_temp, 1.e-12));
      }

      THEN("MinimumPressure is unchanged by shift/scale") {
        REQUIRE(isClose(eos.MinimumPressure(), base_min_pres, 1.e-12));
      }

      THEN("MaximumPressureAtTemperature is unchanged by shift/scale") {
        REQUIRE(isClose(eos.MaximumPressureAtTemperature(0.0), base_max_pres, 1.e-12));
      }

      THEN("RhoPmin returns the scaled minimum density") {
        REQUIRE(isClose(eos.RhoPmin(0.0), base_rho_pmin / scale, 1.e-12));
      }
    }

    AND_GIVEN("A UnitSystem") {
      constexpr Real rho_unit = 2.0;
      constexpr Real sie_unit = 3.0;
      constexpr Real temp_unit = 4.0;
      constexpr Real press_unit = rho_unit * sie_unit;
      auto us =
          UnitSystem<BoundedGas>(BoundedGas(gm1, Cv), rho_unit, sie_unit, temp_unit);

      THEN("The unit system propagates bounds correctly") {
        REQUIRE(isClose(us.MinimumDensity(), base_min_rho / rho_unit, EPS));
        REQUIRE(isClose(us.MinimumTemperature(), base_min_temp / temp_unit, EPS));
        REQUIRE(isClose(us.MaximumDensity(), base_max_rho / rho_unit, EPS));
        REQUIRE(isClose(us.MinimumPressure(), base_min_pres / press_unit, EPS));
        REQUIRE(isClose(us.MaximumPressureAtTemperature(0.0), base_max_pres / press_unit,
                        EPS));
        REQUIRE(isClose(us.RhoPmin(0.0), base_rho_pmin / rho_unit, EPS));
      }
    }
  }
}

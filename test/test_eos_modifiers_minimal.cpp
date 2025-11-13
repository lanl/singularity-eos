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
using singularity::Gruneisen;
using singularity::IdealGas;
using singularity::RelativisticEOS;
using singularity::ScaledEOS;
using singularity::ShiftedEOS;
using singularity::UnitSystem;
using singularity::Variant;

SCENARIO("Relativistic EOS", "[EOSBuilder][RelativisticEOS][IdealGas]") {
  using EOS = Variant<IdealGas, RelativisticEOS<IdealGas>>;
  GIVEN("Parameters for an ideal gas") {
    constexpr Real Cv = 2.0;
    constexpr Real gm1 = 0.5;
    WHEN("We construct a relativistic IdealGas with EOSBuilder") {
      constexpr Real cl = 1;
      EOS eos = IdealGas(gm1, Cv);
      eos = Modify<RelativisticEOS>(eos, cl);
      THEN("The EOS has finite sound speeds") {
        constexpr Real rho = 1e3;
        constexpr Real sie = 1e3;
        Real bmod = eos.BulkModulusFromDensityInternalEnergy(rho, sie);
        Real cs2 = bmod / rho;
        REQUIRE(cs2 < 1);
      }
    }
  }
}

SCENARIO("EOS Unit System", "[EOSBuilder][UnitSystem][IdealGas][Gruneisen]") {
  using EOS = Variant<IdealGas, UnitSystem<IdealGas>, Gruneisen, UnitSystem<Gruneisen>>;
  GIVEN("Parameters for an ideal gas") {
    constexpr Real Cv = 2.0;
    constexpr Real gm1 = 0.5;
    GIVEN("Units with a thermal unit system") {
      constexpr Real rho_unit = 1e1;
      constexpr Real sie_unit = 1e-1;
      constexpr Real temp_unit = 123;
      WHEN("We construct an IdealGas with EOSBuilder") {
        EOS eos = IdealGas(gm1, Cv);
        eos = Modify<UnitSystem>(eos, rho_unit, sie_unit, temp_unit);
        THEN("Units cancel out for an ideal gas") {
          Real rho = 1e3;
          Real sie = 1e3;
          Real P = eos.PressureFromDensityInternalEnergy(rho, sie);
          Real Ptrue = gm1 * rho * sie;
          REQUIRE(std::abs(P - Ptrue) / Ptrue < 1e-3);
        }
      }
    }
    GIVEN("Units with length and time units") {
      constexpr Real time_unit = 456;
      constexpr Real length_unit = 1e2;
      constexpr Real mass_unit = 1e6;
      constexpr Real temp_unit = 789;
      WHEN("We construct an IdealGas with EOSBuilder") {
        EOS eos = IdealGas(gm1, Cv);
        eos = Modify<UnitSystem>(eos, eos_units_init::length_time_units_init_tag,
                                 time_unit, mass_unit, length_unit, temp_unit);
        THEN("Units cancel out for an ideal gas") {
          constexpr Real rho = 1e3;
          constexpr Real sie = 1e3;
          constexpr Real Ptrue = gm1 * rho * sie;
          Real P = eos.PressureFromDensityInternalEnergy(rho, sie);
          REQUIRE(std::abs(P - Ptrue) / Ptrue < 1e-3);
        }
      }
    }
  }
  GIVEN("A Gruneisen EOS") {
    constexpr Real EPS = 10 * std::numeric_limits<Real>::epsilon();

    constexpr Real C0 = 0.394e6;
    constexpr Real s1 = 1.489;
    constexpr Real s2 = 0.0;
    constexpr Real s3 = 0.0;
    constexpr Real G0 = 2.02;
    constexpr Real b = 0.47;
    constexpr Real rho0 = 8.93;
    constexpr Real T0 = 297.0;
    constexpr Real P0 = 1.0e6;
    constexpr Real Cv = 8.56796e7;
    EOS eos = Gruneisen(C0, s1, s2, s3, G0, b, rho0, T0, P0, Cv);
    AND_GIVEN("A Gruneisen EOS with thermal units") {
      constexpr Real rho_unit = 1;
      constexpr Real sie_unit = 2;
      constexpr Real temp_unit = 3;
      EOS eos_units = eos.Modify<UnitSystem>(rho_unit, sie_unit, temp_unit);
      WHEN("We call FillEos") {
        Real rho_cgs = 33;
        Real sie_cgs = 1e13;
        Real P_units, T_units, cv_units, bmod_units;
        Real P_cgs, T_cgs, cv_cgs, bmod_cgs;
        const unsigned long output = (thermalqs::pressure | thermalqs::temperature |
                                      thermalqs::bulk_modulus | thermalqs::specific_heat);
        Real rho_units = rho_cgs / rho_unit;
        Real sie_units = sie_cgs / sie_unit;
        eos_units.FillEos(rho_units, T_units, sie_units, P_units, cv_units, bmod_units,
                          output);

        eos.FillEos(rho_cgs, T_cgs, sie_cgs, P_cgs, cv_cgs, bmod_cgs, output);

        THEN("All quantities agree after unit conversion") {
          REQUIRE(isClose(rho_units * rho_unit, rho_cgs, EPS));
          REQUIRE(isClose(sie_units * sie_unit, sie_cgs, EPS));
          REQUIRE(isClose(T_units * temp_unit, T_cgs, EPS));
          REQUIRE(isClose(P_units * rho_unit * sie_unit, P_cgs, EPS));
          REQUIRE(isClose(cv_units * sie_unit / temp_unit, cv_cgs, EPS));
          REQUIRE(isClose(bmod_units * rho_unit * sie_unit, bmod_cgs, EPS));
        }
      }
      WHEN("We call DensityEnergyFromPressureTemperature") {
        Real P_cgs = 1e6;
        Real T_cgs = 298.15;
        Real P_units = P_cgs / (rho_unit * sie_unit);
        Real T_units = T_cgs / temp_unit;
        Real rho_cgs = 1;
        Real rho_units = 1;
        Real sie_cgs = 1;
        Real sie_units = 1;

        eos.DensityEnergyFromPressureTemperature(
            P_cgs, T_cgs, static_cast<Real *>(nullptr), rho_cgs, sie_cgs);
        eos_units.DensityEnergyFromPressureTemperature(
            P_units, T_units, static_cast<Real *>(nullptr), rho_units, sie_units);
        printf("Rho: %.14e * %.14e  == %.14e, diff: %.14e\n", rho_units, rho_unit,
               rho_cgs, rho_units * rho_unit - rho_cgs);
        REQUIRE(isClose(rho_units * rho_unit, rho_cgs, EPS));
        REQUIRE(isClose(sie_units * sie_unit, sie_cgs, EPS));
      }
    }
  }
}

SCENARIO("Serialization of modified EOSs preserves their properties",
         "[ScaledEOS][IdealGas][Serialization]") {
  using EOS = singularity::Variant<ScaledEOS<IdealGas>, IdealGas>;
  GIVEN("A scaled ideal gas object") {
    constexpr Real Cv = 2.0;
    constexpr Real gm1 = 0.5;
    constexpr Real scale = 2.0;

    constexpr Real rho_test = 1.0;
    constexpr Real sie_test = 1.0;
    constexpr Real temp_trivial = sie_test / (Cv);      // = 1./2.
    constexpr Real temp_test = sie_test / (Cv * scale); // = 1./4.
    constexpr Real EPS = 10 * std::numeric_limits<Real>::epsilon();

    ScaledEOS<IdealGas> eos(IdealGas(gm1, Cv), scale);
    REQUIRE(isClose(eos.TemperatureFromDensityInternalEnergy(rho_test, sie_test),
                    temp_test, EPS));
    EOS eos_scaled = eos;

    EOS eos_trivial = ScaledEOS<IdealGas>(IdealGas(gm1, Cv), 1.0);
    REQUIRE(isClose(eos_trivial.TemperatureFromDensityInternalEnergy(rho_test, sie_test),
                    temp_trivial, EPS));

    THEN("The size of the object is larger than just the ideal gas by itself") {
      REQUIRE(eos.SerializedSizeInBytes() > sizeof(IdealGas));
    }

    WHEN("We serialize the EOS") {
      singularity::VectorSerializer<EOS> serializer({eos_scaled, eos_trivial});
      auto [size, data] = serializer.Serialize();
      REQUIRE(size == serializer.SerializedSizeInBytes());
      REQUIRE(size > 0);
      REQUIRE(data != nullptr);

      THEN("We can de-serialize the EOS") {
        singularity::VectorSerializer<EOS> deserializer;
        deserializer.DeSerialize(data);
        REQUIRE(deserializer.Size() == serializer.Size());

        AND_THEN("The de-serialized EOS still evaluates properly") {
          auto eos_new = deserializer.eos_objects[0];
          REQUIRE(
              isClose(eos_new.TemperatureFromDensityInternalEnergy(rho_test, sie_test),
                      temp_test, EPS));
        }
      }

      free(data);
    }

    eos.Finalize();
  }
}

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
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>
#include <ports-of-call/portable_errors.hpp>
#include <singularity-eos/base/serialization_utils.hpp>
#include <singularity-eos/eos/eos.hpp>
#include <singularity-eos/eos/eos_builder.hpp>
#include <singularity-utils/fast-math/logs.hpp>
#include <singularity-utils/root-finding-1d/root_finding.hpp>

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

namespace {
// A cold energy term exposes density-conversion errors hidden by an ideal gas.
template <bool EnergyPreferred = false>
class ModifierTestEOS : public singularity::EosBase<ModifierTestEOS<EnergyPreferred>> {
 public:
  SG_ADD_BASE_CLASS_USINGS(ModifierTestEOS<EnergyPreferred>)
  PORTABLE_FUNCTION void CheckParams() const {}
  static constexpr unsigned long PreferredInput() {
    return thermalqs::density | (EnergyPreferred ? thermalqs::specific_internal_energy
                                                 : thermalqs::temperature);
  }
  template <typename Lambda = Real *>
  PORTABLE_FUNCTION Real InternalEnergyFromDensityTemperature(Real rho, Real temp,
                                                              Lambda && = nullptr) const {
    return rho * rho + 3 * temp - 5;
  }
  template <typename Lambda = Real *>
  PORTABLE_FUNCTION Real TemperatureFromDensityInternalEnergy(Real rho, Real energy,
                                                              Lambda && = nullptr) const {
    return (energy - rho * rho + 5) / 3;
  }
  template <typename Lambda = Real *>
  PORTABLE_FUNCTION Real MinInternalEnergyFromDensity(Real rho,
                                                      Lambda && = nullptr) const {
    return rho * rho - 5;
  }
  template <typename Lambda = Real *>
  PORTABLE_FUNCTION Real PressureFromDensityTemperature(Real rho, Real temp,
                                                        Lambda && = nullptr) const {
    return rho * (temp + rho);
  }
  template <typename Lambda = Real *>
  PORTABLE_FUNCTION Real SpecificHeatFromDensityTemperature(Real, Real,
                                                            Lambda && = nullptr) const {
    return 3;
  }
  template <typename Lambda = Real *>
  PORTABLE_FUNCTION Real EntropyFromDensityTemperature(Real rho, Real temp,
                                                       Lambda && = nullptr) const {
    return rho + 2 * temp;
  }
  template <typename Lambda = Real *>
  PORTABLE_FUNCTION Real BulkModulusFromDensityTemperature(Real rho, Real temp,
                                                           Lambda && = nullptr) const {
    return rho * (temp + 2 * rho) + rho / 3 * (temp + rho - 2 * rho * rho);
  }
  template <typename Lambda = Real *>
  PORTABLE_FUNCTION Real GruneisenParamFromDensityTemperature(Real, Real,
                                                              Lambda && = nullptr) const {
    return Real(1) / 3;
  }
#define MODIFIER_TEST_FROM_ENERGY(NAME)                                                  \
  template <typename Lambda = Real *>                                                    \
  PORTABLE_FUNCTION Real NAME##FromDensityInternalEnergy(                                \
      Real rho, Real energy, Lambda &&lambda = nullptr) const {                          \
    return NAME##FromDensityTemperature(                                                 \
        rho, TemperatureFromDensityInternalEnergy(rho, energy), lambda);                 \
  }
  MODIFIER_TEST_FROM_ENERGY(Pressure)
  MODIFIER_TEST_FROM_ENERGY(SpecificHeat)
  MODIFIER_TEST_FROM_ENERGY(Entropy)
  MODIFIER_TEST_FROM_ENERGY(BulkModulus)
  MODIFIER_TEST_FROM_ENERGY(GruneisenParam)
#undef MODIFIER_TEST_FROM_ENERGY
  template <typename Lambda = Real *>
  PORTABLE_FUNCTION void InternalEnergyFromDensityPressure(Real rho, Real pressure,
                                                           Real &energy,
                                                           Lambda && = nullptr) const {
    energy = InternalEnergyFromDensityTemperature(rho, pressure / rho - rho);
  }
  template <typename Lambda = Real *>
  PORTABLE_FUNCTION void FillEos(Real &rho, Real &temp, Real &energy, Real &pressure,
                                 Real &cv, Real &bulk, unsigned long output,
                                 Lambda && = nullptr) const {
    if (output & thermalqs::density) {
      rho = (std::sqrt(temp * temp + 4 * pressure) - temp) / 2;
    }
    if (output & thermalqs::temperature)
      temp = TemperatureFromDensityInternalEnergy(rho, energy);
    if (output & thermalqs::specific_internal_energy)
      energy = InternalEnergyFromDensityTemperature(rho, temp);
    if (output & thermalqs::pressure)
      pressure = PressureFromDensityTemperature(rho, temp);
    if (output & thermalqs::specific_heat) cv = 3;
    if (output & thermalqs::bulk_modulus)
      bulk = BulkModulusFromDensityTemperature(rho, temp);
  }
  template <typename Lambda = Real *>
  PORTABLE_FUNCTION void
  ValuesAtReferenceState(Real &rho, Real &temp, Real &energy, Real &pressure, Real &cv,
                         Real &bulk, Real &dpde, Real &dvdt, Lambda && = nullptr) const {
    rho = 2;
    temp = 4;
    energy = 11;
    pressure = 12;
    cv = 3;
    bulk = BulkModulusFromDensityTemperature(rho, temp);
    dpde = Real(2) / 3;
    dvdt = Real(1) / 16;
  }
#define MODIFIER_TEST_VECTOR(NAME)                                                       \
  template <typename Space, typename Lambdas,                                            \
            typename = std::enable_if_t<!variadic_utils::has_int_index_v<Space>>>        \
  void NAME(const Space &, const Real *x, const Real *y, Real *out, Real *, int count,   \
            Lambdas &&, singularity::Transform &&transform = {}) const {                 \
    for (int i = 0; i < count; ++i)                                                      \
      out[i] =                                                                           \
          transform.f.get() * NAME(transform.x.get() * x[i], transform.y.get() * y[i]);  \
  }
  MODIFIER_TEST_VECTOR(InternalEnergyFromDensityTemperature)
  MODIFIER_TEST_VECTOR(TemperatureFromDensityInternalEnergy)
  MODIFIER_TEST_VECTOR(PressureFromDensityTemperature)
  MODIFIER_TEST_VECTOR(PressureFromDensityInternalEnergy)
  MODIFIER_TEST_VECTOR(SpecificHeatFromDensityTemperature)
  MODIFIER_TEST_VECTOR(SpecificHeatFromDensityInternalEnergy)
  MODIFIER_TEST_VECTOR(EntropyFromDensityTemperature)
  MODIFIER_TEST_VECTOR(EntropyFromDensityInternalEnergy)
  MODIFIER_TEST_VECTOR(BulkModulusFromDensityTemperature)
  MODIFIER_TEST_VECTOR(BulkModulusFromDensityInternalEnergy)
  MODIFIER_TEST_VECTOR(GruneisenParamFromDensityTemperature)
  MODIFIER_TEST_VECTOR(GruneisenParamFromDensityInternalEnergy)
#undef MODIFIER_TEST_VECTOR
  static unsigned long scratch_size(std::string method, unsigned int count) {
    return method == "MinInternalEnergyFromDensity" ? count * sizeof(Real) : 0;
  }
  static unsigned long max_scratch_size(unsigned int count) {
    return count * sizeof(Real);
  }
  template <typename Space, typename Lambdas>
  void MinInternalEnergyFromDensity(const Space &, const Real *rhos, Real *out,
                                    Real *scratch, int count, Lambdas &&,
                                    singularity::Transform &&transform = {}) const {
    for (int i = 0; i < count; ++i) {
      scratch[i] = transform.x.get() * rhos[i];
      out[i] = transform.f.get() * MinInternalEnergyFromDensity(scratch[i]);
    }
  }
};

// rho_user=rho_factor*rho_base, e_user=energy_factor*e_base+shift,
// T_user=temp_factor*T_base. Compare against the analytic base at several states.
template <typename EOS>
void checkModifierConversions(EOS &eos, Real rho_factor, Real energy_factor,
                              Real temp_factor, Real shift = 0) {
  const ModifierTestEOS<> base;
  constexpr int count = 3;
  Real rhos[count], temps[count], energies[count], result[count];
  for (int i = 0; i < count; ++i) {
    const Real rho = Real(0.7) + i, temp = 2 + i;
    rhos[i] = rho_factor * rho;
    temps[i] = temp_factor * temp;
    energies[i] =
        energy_factor * base.InternalEnergyFromDensityTemperature(rho, temp) + shift;
    CHECK(isClose(eos.InternalEnergyFromDensityTemperature(rhos[i], temps[i]),
                  energies[i], 1.e-12));
    CHECK(isClose(eos.TemperatureFromDensityInternalEnergy(rhos[i], energies[i]),
                  temps[i], 1.e-12));
    Real r = rhos[i], t = 0, e = energies[i], p = 0, c = 0, b = 0;
    eos.FillEos(r, t, e, p, c, b,
                thermalqs::temperature | thermalqs::pressure | thermalqs::specific_heat |
                    thermalqs::bulk_modulus);
    CHECK(e == energies[i]);
    CHECK(isClose(t, temps[i], 1.e-12));
    CHECK(isClose(c, 3 * energy_factor / temp_factor, 1.e-12));
    CHECK(isClose(
        p, rho_factor * energy_factor * base.PressureFromDensityTemperature(rho, temp),
        1.e-12));
    e = -999;
    eos.FillEos(r, t, e, p, c, b,
                thermalqs::specific_internal_energy | thermalqs::pressure |
                    thermalqs::specific_heat | thermalqs::bulk_modulus);
    CHECK(isClose(e, energies[i], 1.e-12));
    CHECK(isClose(c, 3 * energy_factor / temp_factor, 1.e-12));
    Real solved_energy = energies[i];
    eos.InternalEnergyFromDensityPressure(rhos[i], p, solved_energy);
    CHECK(isClose(solved_energy, energies[i], 1.e-12));
    r = -999;
    e = -999;
    eos.FillEos(r, t, e, p, c, b,
                thermalqs::density | thermalqs::specific_internal_energy |
                    thermalqs::specific_heat | thermalqs::bulk_modulus);
    CHECK(isClose(r, rhos[i], 1.e-12));
    CHECK(isClose(e, energies[i], 1.e-12));
  }
  std::vector<Real> scratch(
      std::max(std::size_t(1),
               static_cast<std::size_t>(EOS::max_scratch_size(count) / sizeof(Real))));
#define CHECK_MODIFIER_VECTOR(NAME, Y, YFACTOR, YOFFSET, FACTOR, OFFSET)                 \
  eos.NAME(PortsOfCall::Exec::Host(), rhos, Y, result, scratch.data(), count,            \
           singularity::NullIndexer{});                                                  \
  for (int i = 0; i < count; ++i) {                                                      \
    const Real rho = rhos[i] / rho_factor;                                               \
    const Real y = (Y[i] - (YOFFSET)) / (YFACTOR);                                       \
    CHECK(isClose(result[i], (FACTOR) * base.NAME(rho, y) + (OFFSET), 1.e-12));          \
  }
  CHECK_MODIFIER_VECTOR(InternalEnergyFromDensityTemperature, temps, temp_factor, 0,
                        energy_factor, shift)
  CHECK_MODIFIER_VECTOR(TemperatureFromDensityInternalEnergy, energies, energy_factor,
                        shift, temp_factor, 0)
  CHECK_MODIFIER_VECTOR(PressureFromDensityInternalEnergy, energies, energy_factor, shift,
                        rho_factor * energy_factor, 0)
  CHECK_MODIFIER_VECTOR(SpecificHeatFromDensityTemperature, temps, temp_factor, 0,
                        energy_factor / temp_factor, 0)
  CHECK_MODIFIER_VECTOR(SpecificHeatFromDensityInternalEnergy, energies, energy_factor,
                        shift, energy_factor / temp_factor, 0)
  CHECK_MODIFIER_VECTOR(EntropyFromDensityInternalEnergy, energies, energy_factor, shift,
                        energy_factor / temp_factor, 0)
  CHECK_MODIFIER_VECTOR(BulkModulusFromDensityInternalEnergy, energies, energy_factor,
                        shift, rho_factor * energy_factor, 0)
  CHECK_MODIFIER_VECTOR(GruneisenParamFromDensityInternalEnergy, energies, energy_factor,
                        shift, 1, 0)
#undef CHECK_MODIFIER_VECTOR
  const auto required =
      EOS::scratch_size("MinInternalEnergyFromDensity", count) / sizeof(Real);
  constexpr Real guard = 123456;
  std::vector<Real> guarded(required + count + 2, guard);
  eos.MinInternalEnergyFromDensity(PortsOfCall::Exec::Host(), rhos, result,
                                   guarded.data() + 1, count, singularity::NullIndexer{});
  for (int i = 0; i < count; ++i)
    CHECK(isClose(
        result[i],
        energy_factor * base.MinInternalEnergyFromDensity(rhos[i] / rho_factor) + shift,
        1.e-12));
  CHECK(guarded.front() == guard);
  for (std::size_t i = required + 1; i < guarded.size(); ++i)
    CHECK(guarded[i] == guard);
  Real r, t, e, p, c, b, dpde, dvdt;
  eos.ValuesAtReferenceState(r, t, e, p, c, b, dpde, dvdt);
  CHECK(isClose(r, 2 * rho_factor, 1.e-12));
  CHECK(isClose(t, 4 * temp_factor, 1.e-12));
  CHECK(isClose(e, 11 * energy_factor + shift, 1.e-12));
  CHECK(isClose(p, 12 * rho_factor * energy_factor, 1.e-12));
  CHECK(isClose(c, 3 * energy_factor / temp_factor, 1.e-12));
  CHECK(isClose(b,
                base.BulkModulusFromDensityTemperature(2, 4) * rho_factor * energy_factor,
                1.e-12));
  CHECK(isClose(dpde, Real(2) / 3 * rho_factor, 1.e-12));
  CHECK(isClose(dvdt, Real(1) / (16 * rho_factor * temp_factor), 1.e-12));
}
} // namespace

TEST_CASE("Modifier conversions compose without tabulated dependencies",
          "[modifiers][Transform]") {
  using Base = ModifierTestEOS<>;
  UnitSystem<Base> units(Base{}, 2, 4, 5);
  checkModifierConversions(units, 0.5, 0.25, 0.2);
  UnitSystem<ShiftedEOS<Base>> units_shift(ShiftedEOS<Base>(Base{}, 10), 2, 4, 5);
  checkModifierConversions(units_shift, 0.5, 0.25, 0.2, 2.5);
  ShiftedEOS<UnitSystem<Base>> shift_units(UnitSystem<Base>(Base{}, 2, 4, 5), 2.5);
  checkModifierConversions(shift_units, 0.5, 0.25, 0.2, 2.5);
  ScaledEOS<Base> scaled(Base{}, 2);
  checkModifierConversions(scaled, 0.5, 2, 1);
  ScaledEOS<ShiftedEOS<Base>> scaled_shift(ShiftedEOS<Base>(Base{}, 10), 2);
  checkModifierConversions(scaled_shift, 0.5, 2, 1, 20);
  ShiftedEOS<ScaledEOS<Base>> shift_scaled(ScaledEOS<Base>(Base{}, 2), 20);
  checkModifierConversions(shift_scaled, 0.5, 2, 1, 20);
  // PreferredInput must not determine which quantities FillEos actually receives.
  using EnergyBase = ModifierTestEOS<true>;
  ShiftedEOS<EnergyBase> shifted_energy(EnergyBase{}, 10);
  checkModifierConversions(shifted_energy, 1, 1, 1, 10);
  ScaledEOS<EnergyBase> scaled_energy(EnergyBase{}, 2);
  checkModifierConversions(scaled_energy, 0.5, 2, 1);
}

TEST_CASE("Shifted FillEos converts the actual energy input", "[ShiftedEOS][IdealGas]") {
  IdealGas base(0.5, 3);
  ShiftedEOS<IdealGas> eos(IdealGas(0.5, 3), 10);
  Real rho = 2, temp = 0, energy = 16, pressure = 0, cv = 0, bulk = 0;
  eos.FillEos(rho, temp, energy, pressure, cv, bulk,
              thermalqs::temperature | thermalqs::pressure | thermalqs::specific_heat |
                  thermalqs::bulk_modulus);
  CHECK(isClose(temp, 2));
  CHECK(energy == 16);
  CHECK(isClose(pressure, base.PressureFromDensityTemperature(rho, temp)));
  CHECK(isClose(cv, 3));
  energy = -999;
  eos.FillEos(rho, temp, energy, pressure, cv, bulk,
              thermalqs::specific_internal_energy | thermalqs::pressure);
  CHECK(isClose(energy, 16));
}

TEST_CASE("Scaling preserves energy and heat capacity", "[ScaledEOS][Gruneisen]") {
  Gruneisen base(0.394e6, 1.489, 0, 0, 2.02, 0.47, 8.93, 297, 1.e6, 8.56796e7);
  constexpr Real scale = 2, rho = 5, temp = 400;
  ScaledEOS<Gruneisen> eos(Gruneisen(base), scale);
  const Real energy =
      scale * base.InternalEnergyFromDensityTemperature(scale * rho, temp);
  const Real cv = scale * base.SpecificHeatFromDensityTemperature(scale * rho, temp);
  CHECK(isClose(eos.InternalEnergyFromDensityTemperature(rho, temp), energy));
  CHECK(isClose(eos.TemperatureFromDensityInternalEnergy(rho, energy), temp));
  CHECK(isClose(eos.SpecificHeatFromDensityTemperature(rho, temp), cv));
  CHECK(isClose(eos.SpecificHeatFromDensityInternalEnergy(rho, energy), cv));
  CHECK(isClose(eos.MinInternalEnergyFromDensity(rho),
                scale * base.MinInternalEnergyFromDensity(scale * rho)));
  Real r = rho, t = 0, e = energy, p = 0, c = 0, b = 0;
  eos.FillEos(r, t, e, p, c, b,
              thermalqs::temperature | thermalqs::pressure | thermalqs::specific_heat |
                  thermalqs::bulk_modulus);
  CHECK(isClose(t, temp));
  CHECK(e == energy);
  CHECK(isClose(c, cv));
}

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
        deserializer.DeSerialize(data.get());
        REQUIRE(deserializer.Size() == serializer.Size());

        AND_THEN("The de-serialized EOS still evaluates properly") {
          auto eos_new = deserializer.eos_objects[0];
          REQUIRE(
              isClose(eos_new.TemperatureFromDensityInternalEnergy(rho_test, sie_test),
                      temp_test, EPS));
        }
      }
    }

    eos.Finalize();
  }
}

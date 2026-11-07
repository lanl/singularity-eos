// © 2026. Triad National Security, LLC. All rights reserved.
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <utility>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <ports-of-call/portable_arrays.hpp>
#include <singularity-eos/eos/eos.hpp>
#include <singularity-utils/indexable_types.hpp>

#ifdef SINGULARITY_USE_SESAME2SPINER
namespace {
using singularity::SesameEOS;
using namespace singularity::thermalqs;

bool close(Real a, Real b, Real tolerance = 1e-10) {
  return std::abs(a - b) <=
         tolerance * std::max(Real(1), std::max(std::abs(a), std::abs(b)));
}

// Analytic surfaces with a density-dependent cold energy, negative energy, and
// zero-coordinate nodes. Quadratic terms distinguish rational from linear
// interpolation and expose missing rho derivatives in the thermodynamic formulas.
sesame2spiner::NativeMaterial nativeMaterial(Real cv = 3, bool entropy = true) {
  const std::vector<Real> r = {0, 0.4, 1.5, 4, 10};
  const std::vector<Real> t = {0, 1, 3, 8, 20};
  std::vector<Real> p, e, s;
  for (Real rho : r) {
    for (Real temp : t) {
      p.push_back(1.2 * rho * temp + rho * rho);
      e.push_back(cv * temp + rho * rho - 5);
      s.push_back(rho + 2 * temp);
    }
  }
  sesame2spiner::NativeMaterial result;
  result.metadata = {1234, 1, 12, 6};
  result.pressure = sesame2spiner::makeNativeTable(r, t, p);
  result.energy = sesame2spiner::makeNativeTable(r, t, e);
  if (entropy) result.entropy = sesame2spiner::makeNativeTable(r, t, s);
  return result;
}

sesame2spiner::NativeMaterial material(Real cv = 3, bool entropy = true) {
  auto result = nativeMaterial(cv, entropy);
  using Grid = sesame2spiner::TemperatureGuessGrid;
  using Regular = Spiner::RegularGrid1D<Real>;
  const Grid rg{{Regular(0, 1.5, 13), Regular(1.5, 10, 35)}};
  const Grid eg{{Regular(result.energy.min(), result.energy.max(), 65)}};
  std::vector<Real> temperatures;
  for (int ir = 0; ir < static_cast<int>(rg.nPoints()); ++ir) {
    for (int ie = 0; ie < static_cast<int>(eg.nPoints()); ++ie) {
      const Real rho = rg.x(ir);
      temperatures.push_back(
          cv == 0 ? Real(0) : std::max((eg.x(ie) - rho * rho + 5) / cv, Real(0)));
    }
  }
  result.temperature_guess =
      sesame2spiner::makeTemperatureGuessTable(rg, eg, temperatures);
  return result;
}

// Compare modifier behavior to an independent native instance in base units.
// rho_user=rho_factor*rho_base, e_user=energy_factor*e_base+shift,
// T_user=temp_factor*T_base. These modifiers imply P_factor=rho_factor*energy_factor.
template <typename EOS>
void checkModifierState(EOS &eos, const SesameEOS &base, Real rho_factor,
                        Real energy_factor, Real temp_factor, Real shift = 0) {
  constexpr Real native_rho = 0.7, native_temp = 2;
  const Real rho = rho_factor * native_rho, temp = temp_factor * native_temp;
  const Real pressure_factor = rho_factor * energy_factor;
  const Real cv_factor = energy_factor / temp_factor;
  const Real energy =
      energy_factor * base.InternalEnergyFromDensityTemperature(native_rho, native_temp) +
      shift;
  const Real pressure =
      pressure_factor * base.PressureFromDensityTemperature(native_rho, native_temp);
  const Real cv =
      cv_factor * base.SpecificHeatFromDensityTemperature(native_rho, native_temp);
  const Real bulk =
      pressure_factor * base.BulkModulusFromDensityTemperature(native_rho, native_temp);
  const Real entropy =
      cv_factor * base.EntropyFromDensityTemperature(native_rho, native_temp);
  const Real minimum =
      energy_factor * base.MinInternalEnergyFromDensity(native_rho) + shift;
  CHECK(close(eos.InternalEnergyFromDensityTemperature(rho, temp), energy));
  CHECK(close(eos.TemperatureFromDensityInternalEnergy(rho, energy), temp));
  CHECK(close(eos.PressureFromDensityTemperature(rho, temp), pressure));
  CHECK(close(eos.PressureFromDensityInternalEnergy(rho, energy), pressure));
  CHECK(close(eos.SpecificHeatFromDensityTemperature(rho, temp), cv));
  CHECK(close(eos.SpecificHeatFromDensityInternalEnergy(rho, energy), cv));
  CHECK(close(eos.BulkModulusFromDensityInternalEnergy(rho, energy), bulk));
  CHECK(close(eos.EntropyFromDensityInternalEnergy(rho, energy), entropy));
  const Real gamma = base.GruneisenParamFromDensityTemperature(native_rho, native_temp);
  CHECK(close(eos.GruneisenParamFromDensityTemperature(rho, temp), gamma));
  CHECK(close(eos.GruneisenParamFromDensityInternalEnergy(rho, energy), gamma));
  CHECK(close(eos.GibbsFreeEnergyFromDensityTemperature(rho, temp),
              energy_factor * base.GibbsFreeEnergyFromDensityTemperature(native_rho,
                                                                         native_temp) +
                  shift));
  CHECK(close(eos.MinInternalEnergyFromDensity(rho), minimum));
  CHECK(close(eos.MinimumInternalEnergy(),
              energy_factor * base.MinimumInternalEnergy() + shift));
  CHECK(close(eos.MaximumInternalEnergy(),
              energy_factor * base.MaximumInternalEnergy() + shift));
  CHECK(close(eos.MaximumDensity(), rho_factor * base.MaximumDensity()));

  Real r = rho, t = -1, e = energy, p = 0, c = 0, b = 0;
  eos.FillEos(r, t, e, p, c, b,
              temperature | singularity::thermalqs::pressure | specific_heat |
                  bulk_modulus);
  CHECK(close(t, temp));
  CHECK(close(e, energy)); // A known input must retain its original units and shift.
  CHECK(close(p, pressure));
  CHECK(close(c, cv));
  CHECK(close(b, bulk));
  t = temp;
  e = -999;
  p = 0;
  eos.FillEos(r, t, e, p, c, b,
              specific_internal_energy | singularity::thermalqs::pressure |
                  specific_heat | bulk_modulus);
  CHECK(close(e, energy));
  CHECK(close(p, pressure));
  CHECK(close(c, cv));

  Real dr, dt, de, dp, dc, db, ddpde, ddvdt;
  Real br, bt, be, bp, bc, bb, bdpde, bdvdt;
  eos.ValuesAtReferenceState(dr, dt, de, dp, dc, db, ddpde, ddvdt);
  base.ValuesAtReferenceState(br, bt, be, bp, bc, bb, bdpde, bdvdt);
  CHECK(close(dr, rho_factor * br));
  CHECK(close(dt, temp_factor * bt));
  CHECK(close(de, energy_factor * be + shift));
  CHECK(close(dp, pressure_factor * bp));
  CHECK(close(dc, cv_factor * bc));
  CHECK(close(db, pressure_factor * bb));
  CHECK(close(ddpde, pressure_factor / energy_factor * bdpde));
  CHECK(close(ddvdt, bdvdt / (rho_factor * temp_factor)));

  Real j1, j2, j3, j4, k1, k2, k3, k4;
  eos.PTDerivativesFromPreferred(rho, energy, pressure, temp,
                                 static_cast<Real *>(nullptr), j1, j2, j3, j4);
  base.PTDerivativesFromPreferred(
      native_rho, base.InternalEnergyFromDensityTemperature(native_rho, native_temp),
      base.PressureFromDensityTemperature(native_rho, native_temp), native_temp,
      static_cast<Real *>(nullptr), k1, k2, k3, k4);
  CHECK(close(j1, energy_factor / pressure_factor * k1));
  CHECK(close(j2, rho_factor / pressure_factor * k2));
  CHECK(close(j3, energy_factor / temp_factor * k3));
  CHECK(close(j4, rho_factor / temp_factor * k4));

  constexpr int count = 2;
  const Real rhos[count] = {rho, rho}, temps[count] = {temp, temp},
             energies[count] = {energy, energy};
  Real result[count] = {};
  std::vector<Real> scratch(
      std::max(std::size_t(1),
               static_cast<std::size_t>(EOS::max_scratch_size(count) / sizeof(Real))),
      0);
#define SG_CHECK_MODIFIER_VECTOR(NAME, X, Y, EXPECTED)                                   \
  eos.NAME(PortsOfCall::Exec::Host(), static_cast<const Real *>(X),                      \
           static_cast<const Real *>(Y), static_cast<Real *>(result), scratch.data(),    \
           count, singularity::NullIndexer{});                                           \
  for (Real value : result)                                                              \
    CHECK(close(value, EXPECTED));
  SG_CHECK_MODIFIER_VECTOR(InternalEnergyFromDensityTemperature, rhos, temps, energy)
  SG_CHECK_MODIFIER_VECTOR(TemperatureFromDensityInternalEnergy, rhos, energies, temp)
  SG_CHECK_MODIFIER_VECTOR(PressureFromDensityTemperature, rhos, temps, pressure)
  SG_CHECK_MODIFIER_VECTOR(PressureFromDensityInternalEnergy, rhos, energies, pressure)
  SG_CHECK_MODIFIER_VECTOR(SpecificHeatFromDensityTemperature, rhos, temps, cv)
  SG_CHECK_MODIFIER_VECTOR(SpecificHeatFromDensityInternalEnergy, rhos, energies, cv)
  SG_CHECK_MODIFIER_VECTOR(BulkModulusFromDensityTemperature, rhos, temps, bulk)
  SG_CHECK_MODIFIER_VECTOR(BulkModulusFromDensityInternalEnergy, rhos, energies, bulk)
  SG_CHECK_MODIFIER_VECTOR(EntropyFromDensityTemperature, rhos, temps, entropy)
  SG_CHECK_MODIFIER_VECTOR(EntropyFromDensityInternalEnergy, rhos, energies, entropy)
  SG_CHECK_MODIFIER_VECTOR(GruneisenParamFromDensityTemperature, rhos, temps, gamma)
  SG_CHECK_MODIFIER_VECTOR(GruneisenParamFromDensityInternalEnergy, rhos, energies, gamma)
#undef SG_CHECK_MODIFIER_VECTOR
  eos.MinInternalEnergyFromDensity(
      PortsOfCall::Exec::Host(), static_cast<const Real *>(rhos),
      static_cast<Real *>(result), scratch.data(), count, singularity::NullIndexer{});
  for (Real value : result)
    CHECK(close(value, minimum));
}
} // namespace

TEST_CASE("SesameEOS composes with unit, shift and scaling modifiers",
          "[SesameEOS][modifiers]") {
  using singularity::ScaledEOS;
  using singularity::ShiftedEOS;
  using singularity::UnitSystem;
  SesameEOS base(material());
  ShiftedEOS<SesameEOS> shifted(SesameEOS(material()), 10);
  checkModifierState(shifted, base, 1, 1, 1, 10);
  shifted.Finalize();
  UnitSystem<SesameEOS> units(SesameEOS(material()), 2, 4, 5);
  checkModifierState(units, base, 0.5, 0.25, 0.2);
  units.Finalize();
  UnitSystem<ShiftedEOS<SesameEOS>> units_shift(
      ShiftedEOS<SesameEOS>(SesameEOS(material()), 10), 2, 4, 5);
  checkModifierState(units_shift, base, 0.5, 0.25, 0.2, 2.5);
  units_shift.Finalize();
  ShiftedEOS<UnitSystem<SesameEOS>> shift_units(
      UnitSystem<SesameEOS>(SesameEOS(material()), 2, 4, 5), 2.5);
  checkModifierState(shift_units, base, 0.5, 0.25, 0.2, 2.5);
  shift_units.Finalize();
  ScaledEOS<SesameEOS> scaled(SesameEOS(material()), 2);
  checkModifierState(scaled, base, 0.5, 2, 1);
  scaled.Finalize();
  ScaledEOS<ShiftedEOS<SesameEOS>> scaled_shift(
      ShiftedEOS<SesameEOS>(SesameEOS(material()), 10), 2);
  checkModifierState(scaled_shift, base, 0.5, 2, 1, 20);
  scaled_shift.Finalize();
  base.Finalize();
}

TEST_CASE("Modified SesameEOS restores its required inverse table",
          "[SesameEOS][modifiers][serialization]") {
  using Modified = singularity::UnitSystem<singularity::ShiftedEOS<SesameEOS>>;
  Modified source(singularity::ShiftedEOS<SesameEOS>(SesameEOS(material()), 10), 2, 4, 5);
  const auto packed = source.Serialize();
  Modified restored;
  restored.DeSerialize(packed.second.get());
  source.Finalize();
  SesameEOS base(material());
  checkModifierState(restored, base, 0.5, 0.25, 0.2, 2.5);
  auto device = restored.GetOnDevice();
  auto *result = static_cast<Real *>(PORTABLE_MALLOC(sizeof(Real)));
  const Real energy = restored.InternalEnergyFromDensityTemperature(0.35, 0.4);
  portableFor(
      "Modified SesameEOS inverse on device", 0, 1, PORTABLE_LAMBDA(int) {
        result[0] = device.TemperatureFromDensityInternalEnergy(0.35, energy);
      });
  Real temperature;
  portableCopyToHost(&temperature, result, sizeof(Real));
  CHECK(close(temperature, 0.4));
  PORTABLE_FREE(result);
  device.Finalize();
  restored.Finalize();
  base.Finalize();
}

TEST_CASE("SesameEOS evaluates rational values and physical derivatives", "[SesameEOS]") {
  SesameEOS eos(material());
  const Real r = 0.7, t = 2;
  const Real expected_p = 1.2 * r * t + r * r;
  const Real expected_e = 3 * t + r * r - 5;
  const auto p = eos.PressureWithDerivatives(r, t);
  const auto e = eos.EnergyWithDerivatives(r, t);
  CHECK(close(p.value, expected_p));
  CHECK(close(e.value, expected_e));
  CHECK(close(p.d_dx2, 1.2 * t + 2 * r));
  CHECK(close(p.d_dx1, 1.2 * r));
  CHECK(close(e.d_dx2, 2 * r));
  CHECK(close(e.d_dx1, 3));
  CHECK(close(eos.SpecificHeatFromDensityTemperature(r, t), 3));
  CHECK(close(eos.GruneisenParamFromDensityTemperature(r, t), 0.4));
  const Real expected_b =
      r * (1.2 * t + 2 * r) + (1.2 * r / 3) * (expected_p / r - r * 2 * r);
  CHECK(close(eos.BulkModulusFromDensityTemperature(r, t), expected_b));
  CHECK(close(eos.EntropyFromDensityTemperature(r, t), r + 2 * t));
  CHECK(eos.MinInternalEnergyFromDensity(r) < 0);
  CHECK(eos.MinimumInternalEnergy() == -5);
  CHECK(eos.MaximumDensity() == 10);
  CHECK(eos.MaximumTemperature() == 20);
  CHECK(eos.matid() == 1234);
  CHECK(eos.MeanAtomicMass() == 12);
  CHECK(eos.MeanAtomicNumber() == 6);
  Real dedp, drdp, dedt, drdt;
  eos.PTDerivativesFromPreferred(r, expected_e, expected_p, t,
                                 static_cast<Real *>(nullptr), dedp, drdp, dedt, drdt);
  CHECK(close(dedp, 2 * r / p.d_dx2));
  CHECK(close(drdp, 1 / p.d_dx2));
  CHECK(close(dedt, 3 - 2 * r * p.d_dx1 / p.d_dx2));
  CHECK(close(drdt, -p.d_dx1 / p.d_dx2));
  const Real h = 1e-5;
  CHECK(close(e.d_dx2, (eos.InternalEnergyFromDensityTemperature(r + h, t) -
                        eos.InternalEnergyFromDensityTemperature(r - h, t)) /
                           (2 * h)));
  eos.Finalize();
}

TEST_CASE("SesameEOS inversions share the forward interpolants", "[SesameEOS][inverse]") {
  SesameEOS eos(material());
  CHECK(eos.TemperatureGuessAvailable());
  // Include a native boundary, interiors, and upper extrapolation. The last
  // case checks secant extrapolation consistency rather than analytic quadratics.
  for (Real t : {Real(0), Real(0.3), Real(2), Real(20), Real(25)}) {
    const Real r = 0.7;
    const Real e = eos.InternalEnergyFromDensityTemperature(r, t);
    const Real p = eos.PressureFromDensityTemperature(r, t);
    CHECK(close(eos.TemperatureFromDensityInternalEnergy(r, e), t));
    CHECK(close(eos.PressureFromDensityInternalEnergy(r, e), p));
    CHECK(close(eos.InternalEnergyFromDensityPressure(r, p), e));
    Real rho_out = 0, e_out = 0;
    eos.DensityEnergyFromPressureTemperature(p, t, static_cast<Real *>(nullptr), rho_out,
                                             e_out);
    CHECK(close(rho_out, r));
    CHECK(close(e_out, e));
    CHECK(close(eos.SpecificHeatFromDensityInternalEnergy(r, e), 3));
    CHECK(close(eos.BulkModulusFromDensityInternalEnergy(r, e),
                eos.BulkModulusFromDensityTemperature(r, t)));
    CHECK(close(eos.GruneisenParamFromDensityInternalEnergy(r, e),
                eos.GruneisenParamFromDensityTemperature(r, t)));
  }
  Real rho = 0.7, temp = 2, e = 0, p = 0, cv = 0, b = 0;
  eos.FillEos(rho, temp, e, p, cv, b,
              specific_internal_energy | pressure | specific_heat | bulk_modulus);
  const Real expected_e = e, expected_p = p;
  temp = -1;
  eos.FillEos(rho, temp, e, p, cv, b,
              temperature | pressure | specific_heat | bulk_modulus);
  CHECK(close(temp, 2));
  CHECK(close(p, expected_p));
  temp = -1;
  e = -1;
  eos.FillEos(rho, temp, e, p, cv, b, temperature | specific_internal_energy);
  CHECK(close(temp, 2));
  CHECK(close(e, expected_e));
  rho = -1;
  eos.FillEos(rho, temp, e, p, cv, b, density | specific_internal_energy);
  CHECK(close(rho, 0.7));
  CHECK(close(e, expected_e));
  eos.Finalize();
}

TEST_CASE("SesameEOS uses table seeds and native gradients in shared Newton",
          "[SesameEOS][inverse][guess]") {
  using Grid = sesame2spiner::TemperatureGuessGrid;
  using Regular = Spiner::RegularGrid1D<Real>;
  const Grid rg{{Regular(0, 10, 3)}}, eg{{Regular(-5, 15, 5)}};
  for (int predictor = 0; predictor < 4; ++predictor) {
    auto data = nativeMaterial(3, true);
    std::vector<Real> temperatures;
    for (int ir = 0; ir < static_cast<int>(rg.nPoints()); ++ir) {
      for (int ie = 0; ie < static_cast<int>(eg.nPoints()); ++ie) {
        const Real e = eg.x(ie);
        // Poor initial values with large, zero, negative, or tiny dT_guess/de.
        // Newton uses the native energy derivative, so these slopes do not
        // control its corrections for queries inside the auxiliary grid.
        temperatures.push_back(predictor == 0   ? 5 + e
                               : predictor == 1 ? 1
                               : predictor == 2 ? 2 - 0.05 * e
                                                : 2 + 1e-6 * e);
      }
    }
    data.temperature_guess =
        sesame2spiner::makeTemperatureGuessTable(rg, eg, temperatures);
    SesameEOS eos(std::move(data));
    const Real rho = 0.7, expected_t = 2;
    const Real energy = eos.InternalEnergyFromDensityTemperature(rho, expected_t);
    Real temperature = -1;
    REQUIRE(eos.TryTemperatureFromDensityInternalEnergy(rho, energy, temperature));
    CHECK(close(temperature, expected_t));
    CHECK(close(eos.InternalEnergyFromDensityTemperature(rho, temperature), energy));
    CHECK(close(eos.SpecificHeatFromDensityInternalEnergy(rho, energy), 3));
    eos.Finalize();
  }
}

TEST_CASE("SesameEOS rejects invalid seed corners without a bracketed fallback",
          "[SesameEOS][inverse][guess]") {
  using Grid = sesame2spiner::TemperatureGuessGrid;
  using Regular = Spiner::RegularGrid1D<Real>;
  const Grid rg{{Regular(0, 10, 2)}}, eg{{Regular(-5, 15, 2)}};
  for (int failure = 0; failure < 4; ++failure) {
    auto data = nativeMaterial(3, false);
    std::vector<Real> temperatures(4, 1);
    temperatures[failure] = -1;
    data.temperature_guess =
        sesame2spiner::makeTemperatureGuessTable(rg, eg, temperatures);
    SesameEOS eos(std::move(data));
    const Real energy = eos.InternalEnergyFromDensityTemperature(0.7, 2);
    Real temperature = -7;
    // A native root exists at T=2, but a seed cell with an invalid corner
    // must fail without bypassing the required auxiliary table.
    CHECK_FALSE(eos.TryTemperatureFromDensityInternalEnergy(0.7, energy, temperature));
    CHECK(temperature == -7);
    using singularity::IndexableTypes::RootStatus;
    singularity::IndexerUtils::VariadicIndexer<RootStatus> status;
    // No valid initial seed exists in this cell, so use minimum temperature.
    CHECK(eos.TemperatureFromDensityInternalEnergy(0.7, energy, status) ==
          eos.MinimumTemperature());
    CHECK(status[RootStatus{}] == static_cast<Real>(RootFinding1D::Status::FAIL));
    eos.Finalize();
  }
}

TEST_CASE("SesameEOS Newton handles zero energy and unreachable targets",
          "[SesameEOS][inverse][Newton]") {
  SesameEOS eos(material());
  constexpr Real rho = 0.7;
  for (Real energy : {Real(-1), Real(-1e-14), Real(0), Real(1e-14), Real(1)}) {
    Real temperature = -7;
    REQUIRE(eos.TryTemperatureFromDensityInternalEnergy(rho, energy, temperature));
    CHECK(close(temperature, (energy - rho * rho + 5) / 3));
    CHECK(std::abs(eos.InternalEnergyFromDensityTemperature(rho, temperature) - energy) <=
          1e-12 * std::max(Real(1), std::abs(energy)));
  }
  Real temperature = -7;
  // This energy requires a negative temperature on the native surface.
  CHECK_FALSE(eos.TryTemperatureFromDensityInternalEnergy(
      rho, eos.MinInternalEnergyFromDensity(rho) - 1, temperature));
  CHECK(temperature == -7);
  CHECK_FALSE(eos.TryTemperatureFromDensityInternalEnergy(
      rho, std::numeric_limits<Real>::quiet_NaN(), temperature));
  CHECK(temperature == -7);
  eos.Finalize();
}

TEST_CASE("SesameEOS returns the seed and reports failed Newton solves like SpinerEOS",
          "[SesameEOS][inverse][RootStatus]") {
  using singularity::IndexableTypes::RootStatus;
  using StatusIndexer = singularity::IndexerUtils::VariadicIndexer<RootStatus>;
  using Grid = sesame2spiner::TemperatureGuessGrid;
  using Regular = Spiner::RegularGrid1D<Real>;
  auto data = nativeMaterial(3, false);
  const Grid rg{{Regular(0, 10, 2)}}, eg{{Regular(-10, 30, 3)}};
  constexpr Real seed = 7, rho = 0.7;
  data.temperature_guess =
      sesame2spiner::makeTemperatureGuessTable(rg, eg, std::vector<Real>(6, seed));
  SesameEOS eos(std::move(data));
  const Real unreachable = eos.MinInternalEnergyFromDensity(rho) - 1;
  const Real reachable = eos.InternalEnergyFromDensityTemperature(rho, 2);
  const Real failed = static_cast<Real>(RootFinding1D::Status::FAIL);
  const Real succeeded = static_cast<Real>(RootFinding1D::Status::SUCCESS);
  StatusIndexer status;
  CHECK(eos.TemperatureFromDensityInternalEnergy(rho, unreachable, status) == seed);
  CHECK(status[RootStatus{}] == failed);
  CHECK(eos.PressureFromDensityInternalEnergy(rho, unreachable, status) ==
        eos.PressureFromDensityTemperature(rho, seed));
  CHECK(status[RootStatus{}] == failed);
  CHECK(close(eos.TemperatureFromDensityInternalEnergy(rho, reachable, status), 2));
  CHECK(status[RootStatus{}] == succeeded); // A later success clears an old failure.
  CHECK(eos.TemperatureFromDensityInternalEnergy(rho, unreachable) == seed);
  Real numeric_lambda[] = {-19};
  CHECK(eos.TemperatureFromDensityInternalEnergy(
            rho, unreachable, static_cast<Real *>(numeric_lambda)) == seed);
  CHECK(numeric_lambda[0] == -19); // RootStatus has no numeric lambda slot.

  const Real rhos[] = {rho, rho}, energies[] = {unreachable, reachable};
  Real temperatures[2] = {};
  StatusIndexer statuses[2];
  eos.TemperatureFromDensityInternalEnergy(
      PortsOfCall::Exec::Host(), static_cast<const Real *>(rhos),
      static_cast<const Real *>(energies), static_cast<Real *>(temperatures),
      static_cast<Real *>(nullptr), 2, static_cast<StatusIndexer *>(statuses));
  CHECK(temperatures[0] == seed);
  CHECK(close(temperatures[1], 2));
  CHECK(statuses[0][RootStatus{}] == failed);
  CHECK(statuses[1][RootStatus{}] == succeeded);

  // The fallback temperature and named status pass through composed modifiers.
  singularity::UnitSystem<singularity::ShiftedEOS<SesameEOS>> modified(
      singularity::ShiftedEOS<SesameEOS>(SesameEOS(eos), 10), 2, 4, 5);
  CHECK(close(modified.TemperatureFromDensityInternalEnergy(
                  rho / 2, (unreachable + 10) / 4, status),
              seed / 5));
  CHECK(status[RootStatus{}] == failed);

  auto device = eos.GetOnDevice();
  auto *result = static_cast<Real *>(PORTABLE_MALLOC(4 * sizeof(Real)));
  portableFor(
      "SesameEOS root status on device", 0, 1, PORTABLE_LAMBDA(int) {
        StatusIndexer lambda;
        result[0] = device.TemperatureFromDensityInternalEnergy(rho, unreachable, lambda);
        result[1] = lambda[RootStatus{}];
        result[2] = device.TemperatureFromDensityInternalEnergy(rho, reachable, lambda);
        result[3] = lambda[RootStatus{}];
      });
  Real host[4];
  portableCopyToHost(host, result, sizeof(host));
  CHECK(host[0] == seed);
  CHECK(host[1] == failed);
  CHECK(close(host[2], 2));
  CHECK(host[3] == succeeded);
  PORTABLE_FREE(result);
  device.Finalize();
  eos.Finalize();
}

TEST_CASE("SesameEOS handles zero density and nonpositive heat capacity", "[SesameEOS]") {
  for (Real heat_capacity : {Real(0), Real(-3)}) {
    SesameEOS eos(material(heat_capacity, false));
    const Real r = 0.7, t = 2;
    CHECK(eos.SpecificHeatFromDensityTemperature(r, t) == 0);
    CHECK(close(eos.BulkModulusFromDensityTemperature(r, t),
                r * eos.PressureWithDerivatives(r, t).d_dx2));
    CHECK(eos.BulkModulusFromDensityTemperature(0, t) == 0);
    CHECK(eos.GruneisenParamFromDensityTemperature(0, t) == 0);
    CHECK_FALSE(eos.EntropyAvailable());
    if (heat_capacity < 0) {
      const Real e = eos.InternalEnergyFromDensityTemperature(r, t);
      CHECK(close(eos.TemperatureFromDensityInternalEnergy(r, e), t));
    } else {
      // A constant thermal branch uses the temperature table's zero seed.
      const Real e = eos.InternalEnergyFromDensityTemperature(r, t);
      CHECK(close(eos.TemperatureFromDensityInternalEnergy(r, e), 0));
    }
    using singularity::IndexableTypes::RootStatus;
    singularity::IndexerUtils::VariadicIndexer<RootStatus> status;
    const Real fallback = eos.TemperatureFromDensityInternalEnergy(r, 1e50, status);
    CHECK(std::isfinite(fallback));
    CHECK(fallback >= 0);
    CHECK(status[RootStatus{}] == static_cast<Real>(RootFinding1D::Status::FAIL));
    const Real reachable = eos.InternalEnergyFromDensityTemperature(r, t);
    eos.TemperatureFromDensityInternalEnergy(r, reachable, status);
    CHECK(status[RootStatus{}] == static_cast<Real>(RootFinding1D::Status::SUCCESS));
#ifdef PORTABILITY_STRATEGY_NONE
    CHECK_THROWS(eos.EntropyFromDensityTemperature(r, t));
#endif
    eos.Finalize();
  }
}

TEST_CASE("SesameEOS selects the branch seeded by the required temperature table",
          "[SesameEOS][inverse]") {
  const std::vector<Real> r = {0, 0.4, 2, 7};
  const std::vector<Real> t = {0, 1, 3, 5, 7, 9};
  std::vector<Real> p, e;
  for (Real rho : r) {
    for (Real temp : t) {
      p.push_back(rho * (1 + temp));
      e.push_back((temp - 4) * (temp - 4));
    }
  }
  using Grid = sesame2spiner::TemperatureGuessGrid;
  using Regular = Spiner::RegularGrid1D<Real>;
  const Grid rg{{Regular(0, 7, 3)}}, eg{{Regular(1, 25, 7)}};
  for (int branch : {-1, 0, 1}) {
    sesame2spiner::NativeMaterial seeded;
    seeded.metadata = {1234, 1, 12, 6};
    seeded.pressure = sesame2spiner::makeNativeTable(r, t, p);
    seeded.energy = sesame2spiner::makeNativeTable(r, t, e);
    std::vector<Real> temperatures;
    for (int ir = 0; ir < static_cast<int>(rg.nPoints()); ++ir)
      for (int ie = 0; ie < static_cast<int>(eg.nPoints()); ++ie)
        temperatures.push_back(std::max(Real(0), 4 + branch * std::sqrt(eg.x(ie))));
    seeded.temperature_guess =
        sesame2spiner::makeTemperatureGuessTable(rg, eg, temperatures);
    SesameEOS eos(std::move(seeded));
    Real temperature = -7;
    // e=6.25 is between auxiliary nodes, so the branch seeds are approximate.
    // The native quadratic derivative must be updated during Newton iteration.
    if (branch != 0) {
      REQUIRE(eos.TryTemperatureFromDensityInternalEnergy(0.4, 6.25, temperature));
      CHECK(close(temperature, 4 + branch * Real(2.5)));
      CHECK(close(eos.InternalEnergyFromDensityTemperature(0.4, temperature), 6.25));
    } else {
      // A seed at the energy minimum has zero native derivative. Roots exist,
      // but Newton must fail instead of searching for a bracket or another seed.
      CHECK_FALSE(eos.TryTemperatureFromDensityInternalEnergy(0.4, 6.25, temperature));
      CHECK(temperature == -7);
      using singularity::IndexableTypes::RootStatus;
      singularity::IndexerUtils::VariadicIndexer<RootStatus> status;
      CHECK(eos.TemperatureFromDensityInternalEnergy(0.4, 6.25, status) == 4);
      CHECK(status[RootStatus{}] == static_cast<Real>(RootFinding1D::Status::FAIL));
    }
    eos.Finalize();
  }
}

TEST_CASE("SesameEOS restores table values and grid caches",
          "[SesameEOS][serialization]") {
  for (bool with_entropy : {false, true}) {
    SesameEOS eos(material(3, with_entropy));
    const auto packed = eos.Serialize();
    const std::size_t bytes = eos.DynamicMemorySizeInBytes();
    SesameEOS restored;
    CHECK(restored.DeSerialize(packed.second.get()) == packed.first);
    std::vector<char> shared(bytes);
    SesameEOS shared_eos;
    CHECK(shared_eos.DeSerialize(packed.second.get(), {shared.data(), true}) ==
          packed.first);
    eos.Finalize();
    CHECK(close(restored.PressureFromDensityTemperature(0.7, 2),
                1.2 * 0.7 * 2 + 0.7 * 0.7));
    CHECK(close(shared_eos.SpecificHeatFromDensityTemperature(0.7, 2), 3));
    CHECK(restored.TemperatureGuessAvailable());
    CHECK(shared_eos.TemperatureGuessAvailable());
    const Real energy = restored.InternalEnergyFromDensityTemperature(0.7, 2);
    CHECK(close(restored.TemperatureFromDensityInternalEnergy(0.7, energy), 2));
    CHECK(close(shared_eos.TemperatureFromDensityInternalEnergy(0.7, energy), 2));
    CHECK(restored.EntropyAvailable() == with_entropy);
    CHECK(restored.DynamicMemorySizeInBytes() == bytes);
    restored.Finalize();
    shared_eos.Finalize();
  }
}

TEST_CASE("SesameEOS supports device queries and host vector dispatch",
          "[SesameEOS][device]") {
  SesameEOS eos(material());
  auto device = eos.GetOnDevice();
  auto *device_result = static_cast<Real *>(PORTABLE_MALLOC(4 * sizeof(Real)));
  portableFor(
      "SesameEOS device queries", 0, 1, PORTABLE_LAMBDA(int) {
        const Real e = device.InternalEnergyFromDensityTemperature(0.7, 2);
        device_result[0] = device.PressureFromDensityTemperature(0.7, 2);
        device_result[1] = device.SpecificHeatFromDensityTemperature(0.7, 2);
        device_result[2] = device.TemperatureFromDensityInternalEnergy(0.7, e);
        device_result[3] = device.EntropyFromDensityTemperature(0.7, 2);
      });
  Real host_result[4];
  portableCopyToHost(host_result, device_result, 4 * sizeof(Real));
  CHECK(close(host_result[0], 1.2 * 0.7 * 2 + 0.7 * 0.7));
  CHECK(close(host_result[1], 3));
  CHECK(close(host_result[2], 2));
  CHECK(close(host_result[3], 4.7));
  PORTABLE_FREE(device_result);
  device.Finalize();
  Real r[2] = {0.7, 1}, t[2] = {2, 3}, p[2] = {};
  eos.PressureFromDensityTemperature(PortsOfCall::Exec::Host(), r, t, p, 2,
                                     singularity::NullIndexer{});
  CHECK(close(p[0], eos.PressureFromDensityTemperature(r[0], t[0])));
  CHECK(close(p[1], eos.PressureFromDensityTemperature(r[1], t[1])));
  singularity::Transform transform;
  transform.x.apply(2);
  transform.y.apply(3);
  transform.f.apply(4);
  eos.PressureFromDensityTemperature(PortsOfCall::Exec::Host(), r, t, p,
                                     static_cast<Real *>(nullptr), 2,
                                     singularity::NullIndexer{}, std::move(transform));
  CHECK(close(p[0], 4 * eos.PressureFromDensityTemperature(2 * r[0], 3 * t[0])));
  CHECK(close(p[1], 4 * eos.PressureFromDensityTemperature(2 * r[1], 3 * t[1])));
  singularity::Variant<SesameEOS> variant(eos);
  CHECK(close(variant.SpecificHeatFromDensityTemperature(0.7, 2), 3));
  eos.Finalize();
}

TEST_CASE("SesameEOS loads real material without retaining EOSPAC handles",
          "[SesameEOS][external]") {
  const char *filename = std::getenv("SINGULARITY_SESAME_FILE");
  const char *id = std::getenv("SINGULARITY_SESAME_MATID");
  if (filename == nullptr || id == nullptr) {
    SKIP("Set SINGULARITY_SESAME_FILE and SINGULARITY_SESAME_MATID for real-data tests");
  }
  SesameEOS eos(filename, std::stoi(id));
  const Real r = (eos.MinimumDensity() + eos.MaximumDensity()) / 2;
  const Real t = (eos.MinimumTemperature() + eos.MaximumTemperature()) / 2;
  const Real energy = eos.InternalEnergyFromDensityTemperature(r, t);
  CHECK(std::isfinite(eos.PressureFromDensityTemperature(r, t)));
  const auto packed = eos.Serialize();
  SesameEOS restored;
  restored.DeSerialize(packed.second.get());
  eos.Finalize();
  CHECK(close(restored.InternalEnergyFromDensityTemperature(r, t), energy));
  restored.Finalize();
}
#endif // SINGULARITY_USE_SESAME2SPINER

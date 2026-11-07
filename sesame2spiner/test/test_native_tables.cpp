// © 2026. Triad National Security, LLC. All rights reserved.
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <utility>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <eos_Interface.h>
#include <eospac-wrapper/eospac_wrapper.hpp>
#include <sesame2spiner/native_tables.hpp>

namespace {
// Independent EOSPAC reference, always selecting the same explicit data file.
struct ReferenceTable {
  EOS_INTEGER handle = -1;
  void load(const char *filename, int matid, EOS_INTEGER type) {
    EOS_INTEGER count = 1, material = matid, error = EOS_OK;
    eos_CreateTables(&count, &type, &material, &handle, &error);
    std::string name(filename);
    eos_SetDataFileName(&handle, &material, &type, name.data(), &error);
    REQUIRE(EospacWrapper::eosErrorCodesEqual(error, EOS_OK));
    eos_LoadTables(&count, &handle, &error);
    REQUIRE(EospacWrapper::eosErrorCodesEqual(error, EOS_OK));
  }
  ~ReferenceTable() {
    EOS_INTEGER count = 1, error = EOS_OK;
    eos_DestroyTables(&count, &handle, &error);
  }
  std::vector<EOS_REAL> info(EOS_INTEGER flag, std::size_t size) {
    EOS_INTEGER count = 1, error = EOS_OK;
    std::vector<EOS_REAL> result(size);
    eos_GetTableInfo(&handle, &count, &flag, result.data(), &error);
    REQUIRE(EospacWrapper::eosErrorCodesEqual(error, EOS_OK));
    return result;
  }
};
} // namespace

TEST_CASE("Native tables preserve coordinates, nodal values and orientation",
          "[sesame2spiner][native]") {
  const std::vector<Real> r = {0, 0.3, 1.5, 7, 13};
  const std::vector<Real> t = {0, 20, 300, 1000};
  std::vector<Real> values;
  for (Real rho : r)
    for (Real temp : t)
      values.push_back(-3 + 2 * rho + 5 * temp);
  sesame2spiner::NativeMaterial owner;
  owner.pressure = sesame2spiner::makeNativeTable(r, t, values);
  const auto &table = owner.pressure;
  REQUIRE(table.dim(2) == r.size());
  REQUIRE(table.dim(1) == t.size());
  for (int ir = 0; ir < r.size(); ++ir) {
    CHECK(table.range(1).x(ir) == r[ir]);
    for (int it = 0; it < t.size(); ++it) {
      CHECK(table.range(0).x(it) == t[it]);
      CHECK(table(ir, it) == values[ir * t.size() + it]);
    }
  }
  const auto result = table.interpToRealWithGrads(Real(0.7), Real(140));
  CHECK(std::abs(result.value - (-3 + 2 * 0.7 + 5 * 140)) < 1e-10);
  CHECK(std::abs(result.d_dx2 - 2) < 1e-10);
  CHECK(std::abs(result.d_dx1 - 5) < 1e-10);
  auto moved = std::move(owner);
  CHECK(moved.pressure(0, 0) == -3);
}

TEST_CASE("Native table validation rejects malformed input", "[sesame2spiner][native]") {
#ifdef PORTABILITY_STRATEGY_NONE
  const std::vector<Real> axis = {0, 1, 2, 3};
  const std::vector<Real> values(16, -7);
  CHECK_THROWS(sesame2spiner::makeNativeTable({0, 1, 2}, axis, values));
  CHECK_THROWS(sesame2spiner::makeNativeTable({0, 1, 1, 3}, axis, values));
  CHECK_THROWS(sesame2spiner::makeNativeTable({-1, 0, 1, 3}, axis, values));
  CHECK_THROWS(sesame2spiner::makeNativeTable(axis, axis, {1, 2}));
  auto invalid = values;
  invalid[3] = std::numeric_limits<Real>::quiet_NaN();
  CHECK_THROWS(sesame2spiner::makeNativeTable(axis, axis, invalid));
#else
  SKIP("Construction errors abort with this portability backend");
#endif
}

TEST_CASE("Native loader rejects invalid construction arguments",
          "[sesame2spiner][native][errors]") {
#ifdef PORTABILITY_STRATEGY_NONE
  using Catch::Matchers::ContainsSubstring;
  // These checks fail before any EOSPAC call and need no material database.
  CHECK_THROWS_WITH(sesame2spiner::loadNativeMaterial("", 1234),
                    ContainsSubstring("An explicit SESAME filename is required"));
  sesame2spiner::NativeLoadOptions options;
  options.split = static_cast<singularity::TableSplit>(-1);
  CHECK_THROWS_WITH(
      sesame2spiner::loadNativeMaterial("unused-sesame-file", 1234, options),
      ContainsSubstring("Invalid native SESAME table split"));
  CHECK_THROWS_WITH(
      sesame2spiner::loadNativeMaterial("unused-sesame-file", 1234, options),
      ContainsSubstring("unused-sesame-file: material 1234"));
#else
  SKIP("Construction errors abort with this portability backend");
#endif
}

TEST_CASE("Auxiliary temperature tables support piecewise grids and inverse slopes",
          "[sesame2spiner][guess]") {
  using Grid = sesame2spiner::TemperatureGuessGrid;
  using Regular = Spiner::RegularGrid1D<Real>;
  const Grid rg{{Regular(0, 1, 5), Regular(1, 10, 7)}};
  const Grid eg{{Regular(-5, 2, 4), Regular(2, 20, 7)}};
  std::vector<Real> values;
  for (int ir = 0; ir < static_cast<int>(rg.nPoints()); ++ir)
    for (int ie = 0; ie < static_cast<int>(eg.nPoints()); ++ie)
      values.push_back(100 + 2 * rg.x(ir) + 3 * eg.x(ie));
  sesame2spiner::NativeMaterial owner;
  owner.temperature_guess = sesame2spiner::makeTemperatureGuessTable(rg, eg, values);
  const auto &table = owner.temperature_guess;
  CHECK(table.range(1).nGrids() == 2);
  CHECK(table.range(0).nGrids() == 2);
  const auto guess = table.interpToRealWithGrads(Real(2), Real(3));
  CHECK(std::abs(guess.value - 113) < 1e-12);
  CHECK(std::abs(guess.d_dx2 - 2) < 1e-12);
  CHECK(std::abs(guess.d_dx1 - 3) < 1e-12);
  auto moved = std::move(owner);
  CHECK(moved.temperature_guess(0, 0) == values[0]);
#ifdef PORTABILITY_STRATEGY_NONE
  CHECK_THROWS(sesame2spiner::makeTemperatureGuessTable(Grid{}, eg, values));
  CHECK_THROWS(sesame2spiner::makeTemperatureGuessTable(rg, eg, {}));
  values[0] = -2;
  CHECK_THROWS(sesame2spiner::makeTemperatureGuessTable(rg, eg, values));
#endif
}

TEST_CASE("Native loader reads an explicitly selected SESAME material",
          "[sesame2spiner][native][external]") {
  const char *filename = std::getenv("SINGULARITY_SESAME_FILE");
  const char *id = std::getenv("SINGULARITY_SESAME_MATID");
  if (filename == nullptr || id == nullptr) {
    SKIP("Set SINGULARITY_SESAME_FILE and SINGULARITY_SESAME_MATID for real-data tests");
  }
  const int matid = std::stoi(id);
  auto material = sesame2spiner::loadNativeMaterial(filename, matid);
  REQUIRE(material.metadata.matid == matid);
  REQUIRE(material.metadata.meanAtomicMass > 0);
  REQUIRE(material.pressure.rank() == 2);
  REQUIRE(material.energy.rank() == 2);
  REQUIRE(material.temperature_guess.rank() == 2);
  CHECK(material.temperature_guess.dim(2) == 2 * (material.energy.dim(2) - 1) + 1);
  CHECK(material.temperature_guess.dim(1) == 2 * (material.energy.dim(1) - 1) + 1);
  const auto &energy = material.energy;
  for (int ir = 0; ir < energy.dim(2); ++ir) {
    for (int it = 0; it < energy.dim(1); ++it) {
      const Real expected = energy(ir, it);
      const auto result =
          energy.interpToRealWithGrads(energy.range(1).x(ir), energy.range(0).x(it));
      CHECK(std::abs(result.value - expected) <=
            1e-12 * std::max(Real(1), std::abs(expected)));
    }
  }
#ifdef PORTABILITY_STRATEGY_NONE
  CHECK_THROWS(sesame2spiner::loadNativeMaterial(filename, -1234567));
  CHECK_THROWS(sesame2spiner::loadNativeMaterial("/nonexistent/sesame-material", matid));
#endif

  // Compare rational values and BOTH gradients, with independent output-unit
  // conversion. Tolerances are provisional until run on production material data.
  const EOS_INTEGER types[2] = {EOS_Pt_DT, EOS_Ut_DT};
  const sesame2spiner::NativeTable *tables[2] = {&material.pressure, &material.energy};
  for (int field = 0; field < 2; ++field) {
    ReferenceTable reference;
    reference.load(filename, matid, types[field]);
    const auto &table = *tables[field];
    REQUIRE(reference.info(EOS_NR, 1)[0] == table.dim(2));
    REQUIRE(reference.info(EOS_NT, 1)[0] == table.dim(1));
    const auto native_r = reference.info(EOS_R_Array, table.dim(2));
    const auto native_t = reference.info(EOS_T_Array, table.dim(1));
    const auto native_f = reference.info(EOS_F_Array, table.size());
    const Real output_scale = field == 0 ? EospacWrapper::pressureFromSesame(1)
                                         : EospacWrapper::sieFromSesame(1);
    for (int ir = 0; ir < table.dim(2); ++ir) {
      CHECK(table.range(1).x(ir) == native_r[ir]);
      for (int it = 0; it < table.dim(1); ++it) {
        CHECK(table.range(0).x(it) == native_t[it]);
        CHECK(table(ir, it) == output_scale * native_f[it * table.dim(2) + ir]);
      }
    }
    for (int ir : {0, (table.dim(2) - 2) / 2, table.dim(2) - 2}) {
      for (int it : {0, (table.dim(1) - 2) / 2, table.dim(1) - 2}) {
        for (Real fraction : {Real(0), Real(0.37), Real(1)}) {
          EOS_REAL rho = table.range(1).x(ir) +
                         fraction * (table.range(1).x(ir + 1) - table.range(1).x(ir));
          EOS_REAL temp = table.range(0).x(it) +
                          fraction * (table.range(0).x(it + 1) - table.range(0).x(it));
          EOS_REAL value, drho, dtemp;
          EOS_INTEGER count = 1, error = EOS_OK;
          eos_Interpolate(&reference.handle, &count, &rho, &temp, &value, &drho, &dtemp,
                          &error);
          REQUIRE(EospacWrapper::eosErrorCodesEqual(error, EOS_OK));
          const auto actual = table.interpToRealWithGrads(Real(rho), Real(temp));
          const Real scale = field == 0 ? EospacWrapper::pressureFromSesame(1)
                                        : EospacWrapper::sieFromSesame(1);
          const auto close = [scale](Real a, Real b) {
            b *= scale;
            return std::abs(a - b) <=
                   1e-8 * std::max(Real(1), std::max(std::abs(a), std::abs(b)));
          };
          INFO("field=" << field << " rho=" << rho << " T=" << temp);
          CHECK(close(actual.value, value));
          CHECK(close(actual.d_dx2, drho));
          CHECK(close(actual.d_dx1, dtemp));
        }
      }
    }
  }
  {
    ReferenceTable inverse;
    inverse.load(filename, matid, EOS_T_DUt);
    const auto &guess = material.temperature_guess;
    for (int ir : {0, guess.dim(2) / 2, guess.dim(2) - 1}) {
      for (int ie : {0, guess.dim(1) / 2, guess.dim(1) - 1}) {
        EOS_REAL rho = EospacWrapper::densityToSesame(guess.range(1).x(ir));
        EOS_REAL energy = EospacWrapper::sieToSesame(guess.range(0).x(ie));
        EOS_REAL temperature, drho, de;
        EOS_INTEGER count = 1, error = EOS_OK;
        eos_Interpolate(&inverse.handle, &count, &rho, &energy, &temperature, &drho, &de,
                        &error);
        REQUIRE((EospacWrapper::eosErrorCodesEqual(error, EOS_OK) ||
                 EospacWrapper::eosErrorCodesEqual(error, EOS_INTERP_EXTRAPOLATED)));
        const Real expected = EospacWrapper::temperatureFromSesame(temperature);
        if (std::isfinite(expected))
          CHECK(std::abs(guess(ir, ie) - std::max(expected, Real(0))) <=
                1e-12 * std::max(std::abs(expected), Real(1)));
        else
          CHECK(guess(ir, ie) == -1);
      }
    }
    // Debug diagnostics call CheckExtrap internally. The checked wrapper must
    // still report the original interpolation status, including extrapolation.
    EOS_REAL rho = EospacWrapper::densityToSesame(guess.range(1).min() - 1);
    EOS_REAL energy = EospacWrapper::sieToSesame(guess.range(0).max());
    EOS_REAL temperature, drho, de;
    EOS_INTEGER count = 1, expected_status = EOS_OK, captured_status = EOS_OK;
    eos_Interpolate(&inverse.handle, &count, &rho, &energy, &temperature, &drho, &de,
                    &expected_status);
    EospacWrapper::eosSafeInterpolate(&inverse.handle, count, &rho, &energy, &temperature,
                                      &drho, &de, "temperature-guess status",
                                      EospacWrapper::Verbosity::Debug, captured_status);
    CHECK(EospacWrapper::eosErrorCodesEqual(captured_status, expected_status));
  }
  {
    sesame2spiner::NativeLoadOptions options;
    options.temperature_guess.rho_points = 7;
    options.temperature_guess.sie_points = 11;
    auto configured = sesame2spiner::loadNativeMaterial(filename, matid, options);
    CHECK(configured.temperature_guess.dim(2) == 7);
    CHECK(configured.temperature_guess.dim(1) == 11);
  }
  const char *other_filename = std::getenv("SINGULARITY_SESAME_OTHER_FILE");
  if (other_filename != nullptr) {
    // Provide a second database with different data for the SAME material ID.
    // Loading it after the first database catches EOSPAC file-cache reuse.
    auto other = sesame2spiner::loadNativeMaterial(other_filename, matid);
    bool different = false;
    const sesame2spiner::NativeTable *other_tables[2] = {&other.pressure, &other.energy};
    for (int field = 0; field < 2; ++field) {
      const auto &a = *tables[field];
      const auto &b = *other_tables[field];
      if (a.dim(1) != b.dim(1) || a.dim(2) != b.dim(2))
        different = true;
      else {
        for (int i = 0; i < a.size(); ++i)
          different = different || a(i) != b(i);
        for (int axis = 0; axis < 2; ++axis)
          for (int i = 0; i < static_cast<int>(a.range(axis).nPoints()); ++i)
            different = different || a.range(axis).x(i) != b.range(axis).x(i);
      }
    }
    CHECK(different);
  }
}

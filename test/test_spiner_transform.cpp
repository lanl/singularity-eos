#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream> // debug

#include <ports-of-call/portability.hpp>
#include <singularity-eos/eos/eos_spiner_common.hpp>
#include <singularity-eos/eos/eos_spiner_sie_transforms.hpp>

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#endif

#include <test/eos_unit_test_helpers.hpp>

#ifdef SPINER_USE_HDF5
#include <singularity-eos/base/spiner_table_utils.hpp>
#include <spiner/databox.hpp>
#include <spiner/interpolation.hpp>
#endif

#ifdef SINGULARITY_USE_SPINER

using namespace singularity;
using spiner_common::DataBox;
using spiner_common::Grid_t;
using spiner_common::to_log;

struct TestDataContainer {
  static constexpr Real lRhoOffset = 0.5;
  static constexpr Real lRho_low = 1.0;
  static constexpr Real lRho_high = 10.0;
  static constexpr size_t Npts = 100;

  using RG = Spiner::RegularGrid1D<Real>;
  Grid_t grid{{RG(lRho_low, lRho_high, Npts)}};
  DataBox sieCold{Npts};
  TestDataContainer() {
    sieCold.setRange(0, grid);
    for (int i = 0; i < Npts; ++i) {
      const Real lRho = grid.x(i);
      sieCold(i) = e_cold_fun(lRho);
    }
  }
  static constexpr Real e_cold_fun(Real lRho) { return 3.0 * lRho + 5.0; }
};

SCENARIO("NullTransform behave correctly", "[TransformTest]") {
  TestDataContainer data;

  NullTransform<> nullTransform;

  GIVEN("A density with a known linear cold-curve in log-space") {
    Real rho = 10.0;
    Real lRho = to_log(rho, data.lRhoOffset);
    Real e_actual = 42.0;

    THEN("NullTransform is identity throughout") {
      Real null_out = nullTransform.transform(e_actual, rho);
      Real null_back = nullTransform.inverse(null_out, lRho);
      REQUIRE(isClose(null_out, e_actual, 1e-14));
      REQUIRE(isClose(null_back, e_actual, 1e-14));
    }
  }
}

SCENARIO("ShiftTransform behave correctly", "[TransformTest]") {
  TestDataContainer data;

  ShiftTransform<TestDataContainer> shiftTransform(data);

  GIVEN("A density with a known linear cold-curve in log-space") {
    Real rho = 10.0;
    Real lRho = to_log(rho, data.lRhoOffset);
    Real cold_curve_value = data.e_cold_fun(lRho);
    Real e_actual = 42.0;

    THEN("Transform subtracts cold curve correctly and Inverse adds cold curve back "
         "correctly") {
      Real e_transformed = shiftTransform.transform(e_actual, rho);
      REQUIRE(isClose(e_transformed, e_actual - cold_curve_value, 1e-14));
      Real e_inverse = shiftTransform.inverse(e_transformed, lRho);
      REQUIRE(isClose(e_inverse, e_actual, 1e-14));
    }
  }
}

#endif // SINGULARITY_USE_SPINER

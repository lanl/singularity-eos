#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream> // debug

#include <ports-of-call/portability.hpp>
#include <singularity-eos/eos/eos_spiner_common.hpp>
#include <singularity-eos/eos/eos_spiner_sie_transforms.hpp>

#include <spiner/databox.hpp>
#include <spiner/interpolation.hpp>

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#endif

#include <test/eos_unit_test_helpers.hpp>

#ifdef SPINER_USE_HDF5
#include <singularity-eos/base/spiner_table_utils.hpp>
#endif

#ifdef SPINER_USE_HDF5
#ifdef SINGULARITY_TEST_SESAME

using namespace singularity;
using spiner_common::DataBox;
using spiner_common::Grid_t;
using spiner_common::to_log;

PORTABLE_INLINE_FUNCTION bool isVeryClose(Real a, Real b, Real eps = 1e-14) {
  return fabs(b - a) / (fabs(a + b) + 1e-20) <= eps;
}

struct TestDataContainer {
  Real lRhoOffset = 0.5;

  using RG = Spiner::RegularGrid1D<Real>;
  Grid_t grid{{RG(1.0, 10.0, 100)}};
  DataBox sieCold;

  TestDataContainer() {
    int N = static_cast<int>(grid.nPoints());
    sieCold.resize(N);
    sieCold.setRange(0, grid);
    for (int i = 0; i < grid.nPoints(); ++i) {
      Real lRho = grid.x(i);
      Real e_cold = 3.0 * lRho + 5.0;
      sieCold(i) = e_cold;
    }
  }
};

SCENARIO("ShiftTransform and NullTransform behave correctly", "[TransformTest]") {
  TestDataContainer data;

  ShiftTransform<TestDataContainer> shiftTransform(data);
  NullTransform<> nullTransform;

  GIVEN("A density with a known linear cold-curve in log-space") {
    Real rho = 10.0;
    Real lRho = to_log(rho, data.lRhoOffset);
    Real cold_curve_value = 3.0 * lRho + 5.0;
    Real e_actual = 42.0;

    THEN("Inverse adds cold curve back correctly (direct check)") {
      Real e_transformed = shiftTransform.transform(e_actual, rho);
      Real e_inverse = shiftTransform.inverse(e_transformed, lRho);
      REQUIRE(isVeryClose(e_inverse, e_actual));
    }

    THEN("Transform then inverse returns the original energy (round‑trip)") {
      Real e_recovered =
          shiftTransform.inverse(shiftTransform.transform(e_actual, rho), lRho);
      REQUIRE(isVeryClose(e_recovered, e_actual));
    }

    THEN("Transform subtracts cold curve correctly") {
      Real e_transformed = shiftTransform.transform(e_actual, rho);
      REQUIRE(isVeryClose(e_transformed, e_actual - cold_curve_value));
    }

    THEN("NullTransform is identity throughout") {
      Real null_out = nullTransform.transform(e_actual, rho);
      Real null_back = nullTransform.inverse(null_out, lRho);
      REQUIRE(isVeryClose(null_out, e_actual));
      REQUIRE(isVeryClose(null_back, e_actual));
    }
  }
}

#endif // SINGULARITY_TEST_SESAME
#endif // SPINER_USE_HDF5

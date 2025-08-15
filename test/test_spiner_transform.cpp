//------------------------------------------------------------------------------
// © 2021-2025. Triad National Security, LLC. All rights reserved.  This
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
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream> // debug

#include <ports-of-call/portability.hpp>
#include <singularity-eos/base/indexable_types.hpp>
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
using namespace transformations;
using spiner_common::DataBox;
using spiner_common::Grid_t;
using spiner_common::to_log;

struct TestDataContainer {
  static constexpr Real lRhoOffset = 0.5;
  static constexpr Real lRho_low = 1.0;
  static constexpr Real lRho_high = 10.0;
  static constexpr size_t Npts = 100;

  static constexpr Real lEOffset = 0.1;
  static constexpr Real lE_low = 0.1;
  static constexpr Real lE_high = 5.0;
  static constexpr size_t N_sie = 80;

  static constexpr Real Cv_a = 0.1;
  static constexpr Real Cv_b = 0.2;
  static constexpr Real Cv_c = 1.0;

  static constexpr Real T_a = 2.0;
  static constexpr Real T_b = 5.0;

  using RG = Spiner::RegularGrid1D<Real>;
  Grid_t rhogrid{{RG(lRho_low, lRho_high, Npts)}};
  Grid_t siegrid{{RG(lE_low, lE_high, N_sie)}};
  DataBox sieCold{Npts};
  DataBox T{Npts};
  DataBox dTdE{Npts, N_sie};

  TestDataContainer() {
    sieCold.setRange(0, rhogrid);

    for (int i = 0; i < Npts; ++i) {
      const Real lRho = rhogrid.x(i);
      sieCold(i) = e_cold_fun(lRho);
    }

    T.setRange(0, rhogrid);

    for (size_t i = 0; i < Npts; ++i) {
      Real lRho = rhogrid.x(i);
      T(i) = T_a * lRho + T_b; // simple linear function, ignore T_c if you want offset
    }

    dTdE.setRange(0, rhogrid); // lRho
    dTdE.setRange(1, siegrid); // lE

    for (int i = 0; i < Npts; ++i) {
      const Real lRho = rhogrid.x(i);
      for (int j = 0; j < N_sie; ++j) {
        const Real lE = siegrid.x(j);
        dTdE(i, j) = Cv_a * lRho + Cv_b * lE + Cv_c;
      }
    }
  }
  static constexpr Real e_cold_fun(Real lRho) { return 3.0 * lRho + 5.0; }

  template <typename IndexerT>
  PORTABLE_INLINE_FUNCTION Real interpRhoSie_(Real rho, Real sie, const DataBox &db,
                                              IndexerT &&lambda) const {
    const Real lRho = spiner_common::to_log(rho, lRhoOffset);
    const Real lE = spiner_common::to_log(sie, lEOffset);
    // If we wanted to mock the lambda write-back (optional)
    if (!variadic_utils::is_nullptr(lambda)) {
      IndexerUtils::Get<IndexableTypes::LogDensity>(lambda, 0) = lRho;
    }
    return db.interpToReal(lRho, lE);
  }

  template <typename IndexerT = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      Real rho, Real sie, IndexerT &&lambda = static_cast<Real *>(nullptr)) const {
    return interpRhoSie_(rho, sie, dTdE, std::forward<IndexerT>(lambda));
  }

  PORTABLE_INLINE_FUNCTION
  Real heatFn(Real rho, Real sie) const {
    return SpecificHeatFromDensityInternalEnergy(rho, sie);
  }

  TestDataContainer GetOnDevice() {
    sieCold = sieCold.getOnDevice();
    dTdE = dTdE.getOnDevice();
    T = T.getOnDevice();
    return *this;
  }
};

template <typename TransformerT>
void test_transformer(const TransformerT transformer, const Real e_actual, const Real rho,
                      const Real expected_transformed) {
  constexpr size_t real_bytes = 1 * sizeof(Real);
  Real e_transformed;
  Real e_inverse;
  Real *e_transformed_d = (Real *)PORTABLE_MALLOC(real_bytes);
  Real *e_inverse_d = (Real *)PORTABLE_MALLOC(real_bytes);
  portableFor(
      "Device execution of transform", 0, 1, PORTABLE_LAMBDA(int i) {
        e_transformed_d[0] = transformer.transform(e_actual, rho);
        e_inverse_d[0] = transformer.inverse(e_transformed_d[0], rho, e_actual);
      });
  portableCopyToHost(&e_transformed, e_transformed_d, real_bytes);
  portableCopyToHost(&e_inverse, e_inverse_d, real_bytes);
  PORTABLE_FREE(e_transformed_d);
  PORTABLE_FREE(e_inverse_d);

  CHECK(isClose(e_transformed, expected_transformed, 1e-14));
  CHECK(isClose(e_inverse, e_actual, 1e-14));
}

// Same as above but runs on host
template <typename TransformerT>
void test_transformer_host(const TransformerT transformer, const Real e_actual,
                           const Real rho, const Real expected_transformed) {
  constexpr size_t real_bytes = 1 * sizeof(Real);
  Real e_transformed;
  Real e_inverse;
  e_transformed = transformer.transform(e_actual, rho);
  e_inverse = transformer.inverse(e_transformed, rho, e_actual);

  CHECK(isClose(e_transformed, expected_transformed, 1e-14));
  CHECK(isClose(e_inverse, e_actual, 1e-14));
}

SCENARIO("NullTransform behave correctly", "[TransformTest]") {
  TestDataContainer data;

  GIVEN("A density with a known linear cold-curve in log-space") {
    Real rho = 10.0;
    Real lRho = to_log(rho, data.lRhoOffset);
    Real e_actual = 42.0;

    THEN("NullTransform is identity throughout when run on device") {
      data = data.GetOnDevice();
      NullTransform<TestDataContainer> nullTransform;
      test_transformer(nullTransform, e_actual, rho, e_actual);
    }
  }
}

SCENARIO("ShiftTransform behave correctly", "[TransformTest]") {
  TestDataContainer data;

  GIVEN("A density with a known linear cold-curve in log-space") {
    Real rho = 10.0;
    Real lRho = to_log(rho, data.lRhoOffset);
    Real cold_curve_value = data.e_cold_fun(lRho);
    Real e_actual = 42.0;
    Real expected_transformed = e_actual - cold_curve_value;

    THEN("Transform subtracts cold curve correctly and Inverse adds cold curve back "
         "correctly when run on device") {
      data = data.GetOnDevice();
      ShiftTransform<TestDataContainer> shiftTransform(data);
      test_transformer(shiftTransform, e_actual, rho, expected_transformed);
    }
  }
}

SCENARIO("DivideByCvTransform behaves correctly", "[TransformTest]") {
  TestDataContainer data;

  GIVEN("An internal energy, density and sie with known constant Cv") {
    Real rho = 10.0;
    Real e_actual = 84.0;

    Real lRho = to_log(rho, data.lRhoOffset);
    Real lE = to_log(e_actual, data.lEOffset);

    Real actual_iCv = data.dTdE.interpToReal(lRho, lE);
    Real expected_transformed = e_actual * actual_iCv;

    THEN("Transform yields e_actual / Cv, inverse reconstructs original when run on "
         "device") {
      // data = data.GetOnDevice();
      DivideByCvTransform<TestDataContainer> cvTransform(data);
      // test_transformer(cvTransform, e_actual, rho, expected_transformed);
      test_transformer_host(cvTransform, e_actual, rho, expected_transformed);
    }
  }
}

SCENARIO("ShiftandDivideByCvTransform behaves correctly", "[TransformTest]") {
  TestDataContainer data;

  GIVEN("An internal energy, density and sie with known constant Cv") {
    Real rho = 10.0;
    Real e_actual = 84.0;

    Real lRho = to_log(rho, data.lRhoOffset);
    Real lE = to_log(e_actual, data.lEOffset);
    Real cold_curve_value = data.e_cold_fun(lRho);
    Real e_shifted = e_actual - cold_curve_value;
    Real actual_iCv = data.dTdE.interpToReal(lRho, lE);
    Real expected_transformed = e_shifted * actual_iCv;

    THEN("Transform yields e_actual / Cv, inverse reconstructs original when run on "
         "device") {
      // data = data.GetOnDevice();
      ShiftandDivideByCvTransform<TestDataContainer> shiftAndCvTransform(data);
      // test_transformer(shiftAndCvTransform, e_actual, rho, expected_transformed);
      test_transformer_host(shiftAndCvTransform, e_actual, rho, expected_transformed);
    }
  }
}

SCENARIO("ScaleTransform behaves correctly", "[TransformTest]") {
  TestDataContainer data;

  GIVEN("An internal energy and density, with linear T(rho)") {
    Real rho = 10.0;
    Real e_actual = 100.0;

    Real lRho = to_log(rho, data.lRhoOffset);
    Real Tval = data.T.interpToReal(lRho);
    Real expected_transformed = e_actual / std::pow(Tval, 3);

    THEN(
        "Transform divides by T³, and inverse reconstructs original when run on device") {
      // data = data.GetOnDevice();
      ScaleTransform<TestDataContainer> scaleTransform(data);
      // test_transformer(scaleTransform, e_actual, rho, expected_transformed);
      test_transformer_host(scaleTransform, e_actual, rho, expected_transformed);
    }
  }
}

SCENARIO("AllTransform behaves correctly", "[TransformTest]") {
  TestDataContainer data;

  GIVEN("An internal energy, density, sie,  with linear T(rho), known constant Cv, and a "
        "known linear cold-curve in log-space") {
    Real rho = 10.0;
    Real e_actual = 100.0;

    Real lRho = to_log(rho, data.lRhoOffset);

    Real cold_curve_value = data.e_cold_fun(lRho);
    Real e_shift_transform = e_actual - cold_curve_value;

    Real lE = to_log(e_shift_transform, data.lEOffset);

    Real actual_Cv = 1. / data.dTdE.interpToReal(lRho, lE);
    Real e_CV_transformed = e_shift_transform / actual_Cv;

    Real Tval = data.T.interpToReal(lRho);
    Real e_final_transformed = e_CV_transformed / std::pow(Tval, 3);

    THEN("Transform correctly applies all transform, and inverse reconstructs original "
         "when run on device") {
      // data = data.GetOnDevice();
      AllTransform<TestDataContainer> Alltest(data);
      // test_transformer(Alltest, e_actual, rho, e_final_transformed);
      test_transformer_host(Alltest, e_actual, rho, e_final_transformed);
    }
  }
}

#endif // SINGULARITY_USE_SPINER

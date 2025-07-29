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
  DataBox sieCold{Npts},T{Npts};
  DataBox cvBox{Npts, N_sie}; 
  

  TestDataContainer() {
    sieCold.setRange(0,rhogrid);


    for (int i = 0; i < Npts; ++i) {
      const Real lRho = rhogrid.x(i);
      sieCold(i) = e_cold_fun(lRho);
    }

    T.setRange(0, rhogrid);

    for (size_t i = 0; i < Npts; ++i) {
      Real lRho = rhogrid.x(i);
      T(i) = T_a * lRho + T_b;  // simple linear function, ignore T_c if you want offset
    }
 
    cvBox.setRange(0, rhogrid); // lRho
    cvBox.setRange(1, siegrid); // lE

    for (int i = 0; i < Npts; ++i) {
      const Real lRho = rhogrid.x(i);
      for (int j = 0; j < N_sie; ++j) {
        const Real lE = siegrid.x(j);
        cvBox(i, j) = Cv_a * lRho + Cv_b * lE + Cv_c;
     }
   }
 }
  static constexpr Real e_cold_fun(Real lRho) { return 3.0 * lRho + 5.0; }


  template <typename IndexerT>
  PORTABLE_INLINE_FUNCTION Real interpRhoSie_(
      Real rho, Real sie, const DataBox &db, IndexerT &&lambda) const {
    const Real lRho = spiner_common::to_log(rho, lRhoOffset);
    const Real lE   = spiner_common::to_log(sie, lEOffset);
    // If we wanted to mock the lambda write-back (optional)
    if (!variadic_utils::is_nullptr(lambda)) {
      IndexerUtils::Get<IndexableTypes::LogDensity>(lambda, 0) = lRho;
    }
    return db.interpToReal(lRho, lE);
  }


  template <typename IndexerT = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      Real rho, Real sie, IndexerT &&lambda = static_cast<Real *>(nullptr)) const {
    return interpRhoSie_(rho, sie, cvBox, std::forward<IndexerT>(lambda));
  }
  
  PORTABLE_INLINE_FUNCTION
     Real heatFn(Real rho, Real sie) const {
     return SpecificHeatFromDensityInternalEnergy(rho, sie);
  }
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

SCENARIO("DivideByCvTransform behaves correctly", "[TransformTest]") {
  TestDataContainer data;
  DivideByCvTransform<TestDataContainer> cvTransform(data);

  GIVEN("An internal energy, density and sie with known constant Cv") {
    Real rho = 10.0;
    Real sie = 42.0;
    Real e_actual = 84.0;

    Real lRho = to_log(rho, data.lRhoOffset);
    Real lE = to_log(sie, data.lEOffset);

    Real actual_Cv = data.cvBox.interpToReal(lRho, lE);
    Real expected_transformed = e_actual / actual_Cv;

    THEN("Transform yields e_actual / Cv, inverse reconstructs original") {
      Real e_transformed = cvTransform.transform(e_actual, rho, sie);
      CHECK(isClose(e_transformed, expected_transformed, 1e-14));

      Real e_inverse = cvTransform.inverse(e_transformed, rho, sie);
      CHECK(isClose(e_inverse, e_actual, 1e-14));
    }
  }
}


SCENARIO("ScaleTransform behaves correctly", "[TransformTest]") {
  TestDataContainer data;
  ScaleTransform<TestDataContainer> scaleTransform(data);

  GIVEN("An internal energy and density, with linear T(rho)") {
    Real rho = 10.0;
    Real e_actual = 100.0;

    Real lRho = to_log(rho, data.lRhoOffset);
    Real Tval = data.T.interpToReal(lRho);
    Real expected_transformed = e_actual / std::pow(Tval, 3);

    THEN("Transform divides by T³, and inverse reconstructs original") {
      Real e_transformed = scaleTransform.transform(e_actual, rho);
      REQUIRE(isClose(e_transformed, expected_transformed, 1e-14));

      Real e_inverse = scaleTransform.inverse(e_transformed, rho);
      REQUIRE(isClose(e_inverse, e_actual, 1e-14));
    }
  }
}

SCENARIO("AllTransform behaves correctly", "[TransformTest]") {
  TestDataContainer data;
  AllTransform<TestDataContainer> Alltest(data);

  GIVEN("An internal energy, density, sie,  with linear T(rho), known constant Cv, and a known linear cold-curve in log-space") {
    Real rho = 10.0;
    Real sie = 42.0;
    Real e_actual = 100.0;

    Real lRho = to_log(rho, data.lRhoOffset);

    Real cold_curve_value = data.e_cold_fun(lRho);
    Real e_shift_transform = e_actual - cold_curve_value;

    Real lE = to_log(e_shift_transform, data.lEOffset);

    Real actual_Cv = data.cvBox.interpToReal(lRho, lE);
    Real e_CV_transformed = e_shift_transform / actual_Cv;


    Real Tval = data.T.interpToReal(lRho);
    Real e_final_transformed = e_CV_transformed/ std::pow(Tval, 3);


    THEN("Transform correctly applies all transform, and inverse reconstructs original") {
      Real e_transformed = Alltest.transform(e_actual, rho);
      REQUIRE(isClose(e_transformed, e_final_transformed, 1e-14));

      Real e_inverse = Alltest.inverse(e_transformed, rho, e_actual);
      REQUIRE(isClose(e_inverse, e_actual, 1e-14));
    }
  }
}




#endif // SINGULARITY_USE_SPINER

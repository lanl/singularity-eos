//------------------------------------------------------------------------------
// © 2021-2026. Triad National Security, LLC. All rights reserved.  This
// program was produced under U.S. Government contract 89233218CNA000001
// for Los Alamos National Laboratory (LANL), which is operated by Triad
// National Security, LLC for the U.S. Department of Energy/National
// Nuclear Security Administration. All rights in the program are
// reserved by Triad National Security, LLC, and the U.S. Department of
// Energy/National Nuclear Security Administration. The Government is
// granted for itself and others acting on its behalf a nonexclusive,
// paid-up, irrevocable worldwide license in this material to reproduce,
// prepare derivative works, distribute copies to the public, perform
// publicly and display publicly, and to permit others to do so.
//------------------------------------------------------------------------------

// This file created with the assistance of generative AI

#include <cmath>
#include <limits>

#include <ports-of-call/portability.hpp>
#include <singularity-eos/base/error_utils.hpp>

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#endif

namespace error_utils = singularity::error_utils;

template <typename T>
void CheckNormalOrZeroClassification(const char *type_name) {
  using limits = std::numeric_limits<T>;
  const T factor = static_cast<T>(error_utils::NORMAL_FACTOR);
  const auto condition = error_utils::is_normal_or_zero{};

  CAPTURE(type_name);

  REQUIRE(condition(T{0}));
  REQUIRE(condition(-T{0}));
  REQUIRE(condition(static_cast<T>(1.382884838243760e+06)));
  REQUIRE(condition(limits::min()));
  REQUIRE(condition(-limits::min()));

  if constexpr (limits::has_denorm != std::denorm_absent) {
    REQUIRE_FALSE(condition(limits::denorm_min()));
    REQUIRE_FALSE(condition(-limits::denorm_min()));
  }

  REQUIRE_FALSE(condition(limits::quiet_NaN()));
  REQUIRE_FALSE(condition(limits::infinity()));
  REQUIRE_FALSE(condition(-limits::infinity()));

  const T safely_bounded = (limits::max() / factor) * T{0.5};
  const T too_large = (limits::max() / factor) * T{2.0};
  REQUIRE(condition(safely_bounded));
  REQUIRE_FALSE(condition(too_large));
}

template <typename T>
void CheckNormalOrZeroClassificationOnDevice(const char *type_name) {
  using limits = std::numeric_limits<T>;
  const T factor = static_cast<T>(error_utils::NORMAL_FACTOR);
  const auto condition = error_utils::is_normal_or_zero{};

  CAPTURE(type_name);

  int nwrong = 0;
  portableReduce(
      "Check error_utils::is_normal_or_zero on device", 0, 1,
      PORTABLE_LAMBDA(const int /*i*/, int &nw) {
        nw += !condition(T{0});
        nw += !condition(-T{0});
        nw += !condition(static_cast<T>(1.382884838243760e+06));
        nw += !condition(limits::min());
        nw += !condition(-limits::min());

        if constexpr (limits::has_denorm != std::denorm_absent) {
          nw += condition(limits::denorm_min());
          nw += condition(-limits::denorm_min());
        }

        nw += condition(limits::quiet_NaN());
        nw += condition(limits::infinity());
        nw += condition(-limits::infinity());
        nw += !condition((limits::max() / factor) * T{0.5});
        nw += condition((limits::max() / factor) * T{2.0});
      },
      nwrong);

  REQUIRE(nwrong == 0);
}

SCENARIO(
    "Error utilities classify normal values, denormals, and overflow-adjacent values",
    "[ErrorUtils]") {
  CheckNormalOrZeroClassification<float>("float");
  CheckNormalOrZeroClassification<double>("double");
}

SCENARIO("Error utilities classify values consistently in device reductions",
         "[ErrorUtils][Device]") {
  CheckNormalOrZeroClassificationOnDevice<float>("float");
  CheckNormalOrZeroClassificationOnDevice<double>("double");
}

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

// This file is made in part with generative AI

#include <cmath>
#include <limits>

#include <singularity-eos/base/robust_utils.hpp>

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#endif

namespace robust = singularity::robust;

SCENARIO("Robust utilities expose scaled floating-point limits", "[RobustUtils]") {
  REQUIRE(robust::SMALL<Real>() == 10 * std::numeric_limits<Real>::min());
  REQUIRE(robust::EPS<Real>() == 10 * std::numeric_limits<Real>::epsilon());
  REQUIRE(robust::min_exp_arg<Real>() ==
          Catch::Approx((std::numeric_limits<Real>::min_exponent - 1) * M_LN2));
  REQUIRE(robust::max_exp_arg<Real>() ==
          Catch::Approx(std::numeric_limits<Real>::max_exponent * M_LN2));
}

SCENARIO("Robust utilities bound scalar values away from unsafe regions",
         "[RobustUtils]") {
  constexpr Real vmin = 0.0;
  constexpr Real vmax = 1.0;

  REQUIRE(robust::make_positive(-1.0) == robust::EPS<Real>());
  REQUIRE(robust::make_positive(2.5) == Catch::Approx(2.5));

  REQUIRE(robust::make_bounded(-1.0, vmin, vmax) == Catch::Approx(vmin + robust::EPS()));
  REQUIRE(robust::make_bounded(0.5, vmin, vmax) == Catch::Approx(0.5));
  REQUIRE(robust::make_bounded(2.0, vmin, vmax) ==
          Catch::Approx(vmax * (1.0 - robust::EPS())));
}

SCENARIO("Robust utilities preserve sign information and regularize ratios",
         "[RobustUtils]") {
  REQUIRE(robust::sgn(-3.0) == -1);
  REQUIRE(robust::sgn(0.0) == 1);
  REQUIRE(robust::sgn(4.0) == 1);
  REQUIRE(robust::sgn(-3) == -1);
  REQUIRE(robust::sgn(0) == 1);
  REQUIRE(robust::sgn(4) == 1);
  REQUIRE(robust::sgn(0u) == 1);
  REQUIRE(robust::sgn(4u) == 1);

  REQUIRE(robust::ratio(12.0, 3.0) == Catch::Approx(4.0));
  REQUIRE(robust::ratio(12.0, -3.0) == Catch::Approx(-4.0));

  const Real ratio_with_zero = robust::ratio(1.0, 0.0);
  REQUIRE(std::isfinite(ratio_with_zero));
  REQUIRE(ratio_with_zero == Catch::Approx(1.0 / robust::SMALL<Real>()));
}

SCENARIO("Robust utilities clamp exponential evaluation to a safe range",
         "[RobustUtils]") {
  REQUIRE(robust::safe_arg_exp(1.5) == Catch::Approx(std::exp(1.5)));
  REQUIRE(robust::safe_arg_exp(robust::min_exp_arg<Real>() - 1.0) == 0.0);
  REQUIRE(std::isinf(robust::safe_arg_exp(robust::max_exp_arg<Real>() + 1.0)));
}

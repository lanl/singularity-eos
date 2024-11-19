//======================================================================
// © 2021-2023. Triad National Security, LLC. All rights reserved.  This
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
//======================================================================

#ifndef _SINGULARITY_EOS_UTILS_FAST_MATH_LOGS_
#define _SINGULARITY_EOS_UTILS_FAST_MATH_LOGS_
#include <cassert>
#include <cmath>
#include <cstdint>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_errors.hpp>
/*
 * These functions are for use when moving to and from the gridding space
 * for Spiner grids.
 *
 * The core idea here is to move onto a grid that's roughly
 * logarithmically spaced. It doesn't actually matter if that grid is
 * EXACTLY logarithmic. Just that it is approximately so.
 * There are however some constraints that do matter:
 *
 * 1. The function we use to translate to "grid space" must be
 * precise, meaning you get the same answer each time.
 *
 * 2. The function must be invertible, meaning you can go back to
 * "linear space." Ideally it is also continuous.
 *
 * 3. The function and its inverse must be fast.
 *
 * To meet these constraints, we split a real number into its mantissa
 * + exponent in base 2. The mantissa is a real number on the interval
 * [0.5, 1) and the exponent is an integer such that
 *
 * mantissa * 2^exponent = number
 *
 * The log in base 2 is then
 *
 * lg(number) = lg(mantissa) + exponent
 *
 * so our goal is to approximate lg(mantissa) satisfying the above
 * criteria. If we do, then we will achieve the desired goals. We have
 * two approaches. The first is a linear approximation of
 * lg(mantissa). The second is a quadratic one. The latter allows for
 * continuity of the derivative.
 *
 * See ArXiv:2206.08957
 */

// TODO(JMM): All the math here is all assuming inputs are
// doubles. Should add single-precision overloads.

namespace singularity {
namespace FastMath {

// TODO(JMM): switch from preprocessor macros to CMake generated C++
namespace Settings {
#ifdef SINGULARITY_USE_TRUE_LOG_GRIDDING
constexpr bool TRUE_LOGS = true;
#else
constexpr bool TRUE_LOGS = false;
#endif // SINGULARITY_USE_TRUE_LOG_GRIDDING
constexpr bool NQT_LOGS = !TRUE_LOGS;

#ifdef SINGULARITY_USE_SINGLE_LOGS
constexpr bool FP32 = true;
#else
constexpr bool FP32 = false;
#endif // SINGULARITY_USE_SINGLE_LOGS
constexpr bool FP64 = !FP32;

#ifdef SINGULARITY_NQT_PORTABLE
constexpr bool NQT_PORTABLE = true;
#else
constexpr bool NQT_PORTABLE = false;
#endif // SINGULARITY_NQT_PORTABLE

#ifdef SINGULARITY_NQT_ORDER_1
constexpr bool NQT_O1 = true;
#else
constexpr bool NQT_O1 = false;
#endif // SINGULARITY_NQT_O1
} // namespace Settings

template <typename T>
PORTABLE_FORCEINLINE_FUNCTION auto sgn(const T &x) {
  return (T(0) < x) - (x < T(0));
}

// TODO(JMM): Add templating and concepts to enforce bit width
// equivalence in input/output types.
PORTABLE_FORCEINLINE_FUNCTION
auto as_int(double f) { return *reinterpret_cast<std::int64_t *>(&f); }
PORTABLE_FORCEINLINE_FUNCTION
auto as_double(std::int64_t i) { return *reinterpret_cast<double *>(&i); }

// Magic numbers constexpr because C++ doesn't constexpr reinterpret casts
// These numbers will be different for different precisions and endianness.
namespace FP64LE { // 64 bit, little endian
constexpr std::int64_t one = 1;
// as_int(1.0) == 2^62 - 2^52
constexpr std::int64_t one_as_int = (one << 62) - (one << 52);
// 1./static_cast<double>(as_int(2.0) - as_int(1.0)) == 2^-52
constexpr double scale_down = 2.22044604925031e-16;
// as_int(2.0) - as_int(1.0) = 2^52, but note the type
constexpr double scale_up = (one << 52);
// 2^52 - 1
constexpr std::int64_t mantissa_mask = (one << 52) - one;
// 2^26 - 1
constexpr std::int64_t low_mask = (one << 26) - 1;
} // namespace FP64LE

// First order interpolation based NQTs
// ----------------------------------------------------------------------
// Reference implementations, however the integer cast implementation
// below is probably faster.
PORTABLE_FORCEINLINE_FUNCTION
double lg_o1_portable(const double x) {
  int e;
  PORTABLE_REQUIRE(x > 0, "log divergent for x <= 0");
  const double m = frexp(x, &e);
  return 2 * (m - 1) + e;
}

PORTABLE_FORCEINLINE_FUNCTION
double pow2_o1_portable(const double x) {
  const int flr = std::floor(x);
  const double remainder = x - flr;
  const double mantissa = 0.5 * (remainder + 1);
  const double exponent = flr + 1;
  return ldexp(mantissa, exponent);
}

// Integer aliased versions
PORTABLE_FORCEINLINE_FUNCTION
double lg_o1_aliased(const double x) {
  using namespace FP64LE;
  PORTABLE_REQUIRE(x > 0, "log divergent for x <= 0");
  return static_cast<double>(as_int(x) - one_as_int) * scale_down;
}

PORTABLE_FORCEINLINE_FUNCTION
double pow2_o1_aliased(const double x) {
  using namespace FP64LE;
  return as_double(static_cast<std::int64_t>(x * scale_up) + one_as_int);
}
// ----------------------------------------------------------------------

// Second-order interpolation based NQTs
// These implementations are due to Peter Hammond
// ----------------------------------------------------------------------
// Portable versions that use frexp/ldexp rather than integer aliasing
PORTABLE_FORCEINLINE_FUNCTION
double lg_o2_portable(const double x) {
  PORTABLE_REQUIRE(x > 0, "log divergent for x <= 0");
  constexpr double four_thirds = 4. / 3.;
  int e;
  const double m = frexp(x, &e);
  return e - four_thirds * (m - 2) * (m - 1);
}

// This version uses the exact formula
PORTABLE_FORCEINLINE_FUNCTION
double pow2_o2_portable(const double x) {
  // log2(mantissa). should go between -1 and 0
  const int flr = std::floor(x);
  const double lm = x - flr - 1;
  const double mantissa = 0.5 * (3 - std::sqrt(1 - 3 * lm));
  const double exponent = flr + 1;
  return ldexp(mantissa, exponent);
}

// Integer aliased/bithacked versions
PORTABLE_FORCEINLINE_FUNCTION
double lg_o2_aliased(const double x) {
  using namespace FP64LE;
  PORTABLE_REQUIRE(x > 0, "log divergent for x <= 0");
  const std::int64_t x_as_int = as_int(x) - one_as_int;
  const std::int64_t frac_as_int = x_as_int & mantissa_mask;
  const std::int64_t frac_high = frac_as_int >> 26;
  const std::int64_t frac_low = frac_as_int & low_mask;
  const std::int64_t frac_squared =
      frac_high * frac_high + ((frac_high * frac_low) >> 25);

  return static_cast<double>(x_as_int + ((frac_as_int - frac_squared) / 3)) * scale_down;
}

PORTABLE_FORCEINLINE_FUNCTION
double pow2_o2_aliased(const double x) {
  using namespace FP64LE;
  constexpr std::int64_t a = 9007199254740992;  // 2 * 2^52
  constexpr double b = 67108864;                // 2^26
  constexpr std::int64_t c = 18014398509481984; // 4 * 2^52
  const std::int64_t x_as_int = static_cast<std::int64_t>(x * scale_up);
  const std::int64_t frac_as_int = x_as_int & mantissa_mask;
  const std::int64_t frac_sqrt =
      static_cast<std::int64_t>(b * std::sqrt(static_cast<double>(c - 3 * frac_as_int)));

  return as_double(x_as_int + a - frac_sqrt - frac_as_int + one_as_int);
}
// ----------------------------------------------------------------------

PORTABLE_FORCEINLINE_FUNCTION
double fastlg(const double x) {
#ifdef SINGULARITY_NQT_ORDER_1
#ifdef SINGULARITY_NQT_PORTABLE
  return lg_o1_portable(x);
#else
  return lg_o1_aliased(x);
#endif // PORTABLE
#else
#ifdef SINGULARITY_NQT_PORTABLE
  return lg_o2_portable(x);
#else
  return lg_o2_aliased(x);
#endif // PORTABLE
#endif // NQT_ORDER
}

PORTABLE_FORCEINLINE_FUNCTION
double fastpow2(const double x) {
#ifdef SINGULARITY_NQT_ORDER_1
#ifdef SINGULARITY_NQT_PORTABLE
  return pow2_o1_portable(x);
#else
  return pow2_o1_aliased(x);
#endif // PORTABLE
#else
#ifdef SINGULARITY_NQT_PORTABLE
  return pow2_o2_portable(x);
#else
  return pow2_o2_aliased(x);
#endif // PORTABLE
#endif // NQT_ORDER
}

PORTABLE_FORCEINLINE_FUNCTION
double lg(const double x) {
  PORTABLE_REQUIRE(x > 0, "log divergent for x <= 0");
  if constexpr (Settings::NQT_LOGS) {
    // Default expression
    return fastlg(x);
  } else {
    if constexpr (Settings::FP64) {
      // double precision is default
      return std::log2(x);
    } else {
      return std::log2f(x);
    }
  }
}

PORTABLE_FORCEINLINE_FUNCTION
double pow2(const double x) {
  if constexpr (Settings::NQT_LOGS) {
    // Default expression
    return fastpow2(x);
  } else {
    if constexpr (Settings::FP64) {
      // double precision is default
      return std::exp2(x);
    } else {
      return std::exp2f(x);
    }
  }
}

PORTABLE_FORCEINLINE_FUNCTION
double log10(const double x) {
  constexpr double LOG2OLOG10 = 0.301029995663981195;
  return LOG2OLOG10 * lg(x);
}

PORTABLE_FORCEINLINE_FUNCTION
double pow10(const double x) {
  constexpr double LOG10OLOG2 = 3.321928094887362626;
  return pow2(LOG10OLOG2 * x);
}

} // namespace FastMath
} // namespace singularity

#endif //  _SINGULARITY_EOS_UTILS_FAST_MATH_LOGS_

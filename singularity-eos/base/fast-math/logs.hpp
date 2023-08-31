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

#include <ports-of-call/portability.hpp>

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
 * To meet these constraints, we approximate a logarithm by using
 * frexp and ldexp to split a real number into its mantissa + exponent
 * in base 2. The mantissa is a real number on the interval [0.5, 1)
 * and the exponent is an integer such that
 *
 * mantissa * 2^exponent = number
 *
 * To truly take the log, one needs to take the log of the mantissa.
 * However, the first term in the logarithm's Taylor series is linear.
 * Thus a linear approximation to log, for small numbers is valid.
 * We use that linear approximation, which turns out to be:
 *
 * 2*(mantissa - 1)
 *
 * So that mantissa = 0.5 returns log_2(0.5) = -1 and mantissa (1)
 * yields 0.
 *
 * Then the approximate log of the whole thing is the approximation of
 * the log of the mantissa, plus the exponent:
 *
 * 2*(mantissa - 1) + exponent
 *
 * Which is a continuous, invertible function that approximates lg
 * only relatively accurately. The absolute difference of this
 * function from lg is about 0.1 at its largest. This translates to a
 * relative difference of at most 25%.
 *
 * However, we don't mind this difference, as we're only interested in
 * generating a grid spacing "close enough" to log spacing. As long as
 * we can go into and out of "grid space" reliably and quickly, we're
 * happy, and this function is an EXACT map into and out of that
 * space.
 *
 * See ArXiv:2206.08957
 */

// TODO(JMM): All the math here is all assuming inputs are
// doubles. Should add single-precision overloads.

namespace singularity {
namespace FastMath {

// Integer Aliased versions. Likely not as portable
// TODO(JMM): Add single-precision versions
PORTABLE_FORCEINLINE_FUNCTION
auto as_long_int(double f) { return *reinterpret_cast<long long int *>(&f); }
PORTABLE_FORCEINLINE_FUNCTION
auto as_double(long long int i) { return *reinterpret_cast<double *>(&i); }

PORTABLE_FORCEINLINE_FUNCTION
double fastlg_aliased(const double x) {
  // Magic numbers constexpr because C++ doesn't constexpr reinterpret casts
  // these are floating point numbers as reinterpreted as integers.
  // as_long_int(1.0)
  constexpr long long int one_as_long_int = 4607182418800017408;
  // 1./static_cast<double>(as_long_int(2.0) - as_long_int(1.0))
  constexpr double scale_down = 2.22044604925031e-16;
  return static_cast<double>(as_long_int(x) - one_as_long_int) * scale_down;
}
PORTABLE_FORCEINLINE_FUNCTION
double fastpow2_aliased(const double x) {
  // Magic numbers constexpr because C++ doesn't constexpr reinterpret casts
  // these are floating point numbers as reinterpreted as integers.
  // as_long_int(1.0)
  constexpr long long int one_as_long_int = 4607182418800017408;
  // as_long_int(2.0) - as_long_int(1.0)
  constexpr double scale_up = 4503599627370496;
  return as_double(static_cast<long long int>(x * scale_up) + one_as_long_int);
}

PORTABLE_FORCEINLINE_FUNCTION
double fastlg(const double x) {
#ifdef SINGULARITY_USE_HIGH_RISK_MATH
  return fastlg_aliased(x);
#else
  int n;
  assert(x > 0 && "log divergent for x <= 0");
#ifndef SINGULARITY_USE_SINGLE_LOGS
  // double precision is default
  const double y = frexp(x, &n);
#else
  // faster but less accurate
  const float y = frexpf((float)x, &n);
#endif // SINGULARITY_USE_SINGLE_LOGS
  return 2 * (y - 1) + n;
#endif // SINGULARITY_USE_HIGH_RISK_MATH
}

PORTABLE_FORCEINLINE_FUNCTION
double fastpow2(const double x) {
#ifdef SINGULARITY_USE_HIGH_RISK_MATH
  return fastpow2_aliased(x);
#else
  const int flr = std::floor(x);
  const double remainder = x - flr;
  const double mantissa = 0.5 * (remainder + 1);
  const double exponent = flr + 1;
#ifndef SINGULARITY_USE_SINGLE_LOGS
  // double precision is default
  return ldexp(mantissa, exponent);
#else
  // Faster but less accurate
  return ldexpf((float)mantissa, exponent);
#endif // SINGULARITY_USE_SINGLE_LOGS
#endif // SINGULARITY_USE_HIGH_RISK_MATH
}

PORTABLE_FORCEINLINE_FUNCTION
double lg(const double x) {

  assert(x > 0 && "log divergent for x <= 0");

#ifndef SINGULARITY_USE_TRUE_LOG_GRIDDING
  // Default expression
  return fastlg(x);
#else
#ifndef SINGULARITY_USE_SINGLE_LOGS
  // double precision is default
  return std::log2(x);
#else
  return std::log2f(x);
#endif // SINGULARITY_USE_SINGLE_LOGS
#endif // SINGULARITY_USE_TRUE_LOG_GRIDDING
}

PORTABLE_FORCEINLINE_FUNCTION
double pow2(const double x) {
#ifndef SINGULARITY_USE_TRUE_LOG_GRIDDING
  // Default expression
  return fastpow2(x);
#else
#ifndef SINGULARITY_USE_SINGLE_LOGS
  // double precision is default
  return std::exp2(x);
#else
  return std::exp2f(x);
#endif // SINGULARITY_USE_SINGLE_LOGS
#endif // SINGULARITY_USE_TRUE_LOG_GRIDDING
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

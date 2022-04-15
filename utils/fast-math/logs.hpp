//======================================================================
// approximate log10, approximately 5x faster than std::log10
// Author: Jonah Miller (jonahm@lanl.gov)
// © 2021-2022. Triad National Security, LLC. All rights reserved.  This
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
#include <cmath>
#include <ports-of-call/portability.hpp>

namespace singularity {
namespace Math {

PORTABLE_FORCEINLINE_FUNCTION
double lg(const double x) {

  int n;
#ifndef SINGULARITY_USE_SINGLE_LOGS
  const double y = frexp(x, &n); // default is double precision
#else
  const float y = frexpf((float)x, &n); // faster but less accurate
#endif // SINGULARITY_USE_SINGLE_LOGS

  return 2*(y - 1) + n;
}

PORTABLE_FORCEINLINE_FUNCTION
double log10(const double x) {
  constexpr double LOG2OLOG10 = 0.301029995663981195;
  return LOG2OLOG10*lg(x);
}

PORTABLE_FORCEINLINE_FUNCTION
double pow2(const double x) {
  const int flr = std::floor(x);
  const double remainder = x - flr;
  const double mantissa = 0.5*(remainder + 1);
  const double exponent = flr + 1;
#ifndef SINGULARITY_USE_SINGLE_LOGS
  return ldexp(mantissa, exponent);
#else
  return ldexpf(mantissa, exponent);
#endif // SINGULARITY_USE_SINGLE_LOGS
}

PORTABLE_FORCEINLINE_FUNCTION
double pow10(const double x) {
  constexpr double LOG10OLOG2 = 3.321928094887362626;
  return pow2(LOG10OLOG2*x);
}

} // namespace Math
} // namespace singularity

#endif //  _SINGULARITY_EOS_UTILS_FAST_MATH_LOGS_

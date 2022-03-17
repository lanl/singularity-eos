//======================================================================
// approximate log10, approximately 5x faster than std::log10
// Author: Jonah Miller (jonahm@lanl.gov)
// © 2021 - 2022. Triad National Security, LLC. All rights reserved.  This
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
double log10(const double x) {

  // const double ILOG10 = 1./std::log(10.0);
  constexpr double ILOG10 = 0.43429448190325176;

  int n{};
  constexpr const double LOG2OLOG10 = 0.301029995663981195;
#ifndef SINGULARITY_USE_SINGLE_LOGS
  const double y = frexp(x, &n); // default is double precision
#else
  const float y = frexpf((float)x, &n); // faster but less accurate
#endif
#ifdef SINGULARITY_FMATH_USE_ORDER_4
  // 4th order approximation to log(x) x in [0.5, 1.0)
  const double expr =
      -2.36697629372863 +
      y * (5.2672858163176 +
           y * (-5.1242620871906 + (2.91607506431871 - 0.69212249971711 * y) * y));
#elif defined(SINGULARITY_FMATH_USE_ORDER_5)
  // 5th order approximation to log(x) x in [0.5, 1.0)
  const double expr =
      -2.5396743497743 +
      y * (6.420250933346 +
           y * (-8.150589926467 +
                y * (6.830564334687 + (-3.192132453902 + 0.6315814621098 * y) * y)));
#else
  // 7the order approximation to log(x) x in [0.5, 1.0)
  const double expr =
      -2.808428971 +
      y * (8.648808309 +
           y * (-15.910426569 +
                y * (21.53943522 +
                     y * (-19.56870772 +
                          y * (11.318698784 + (-3.770730277 + 0.5513512194 * y) * y)))));
#endif // SINGULARITY_FMATH_USE_ORDER
  // log10 x using frexp
  return ILOG10 * expr + n * LOG2OLOG10;
}

} // namespace Math
} // namespace singularity

#endif //  _SINGULARITY_EOS_UTILS_FAST_MATH_LOGS_

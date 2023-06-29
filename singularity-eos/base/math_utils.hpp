//------------------------------------------------------------------------------
// © 2022. Triad National Security, LLC. All rights reserved.  This
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

#ifndef SINGULARITY_EOS_BASE_MATH_UTILS_HPP_
#define SINGULARITY_EOS_BASE_MATH_UTILS_HPP_

#include <ports-of-call/portability.hpp>

namespace singularity {
namespace math_utils {

/*
 * Replaces "SQUARE, CUBE, etc
 * Use ternary operator instead of if constexpr to enable compile-time tail-recursion
 * for various edge cases. Following accelerants are used:
 * x^{2*i} = x^i x^i
 * x^i = i x^{i - 1}
 * x^{-i} = 1/x^i
 */
template <int P>
PORTABLE_FORCEINLINE_FUNCTION constexpr Real pow(const Real x) {
  // defining Ppos is important. Otherwise compiler is not smart
  // enough to terminate recursion.
  constexpr int Ppos = (P < 0) ? -P : P;
  return (P >= 0)
             ? ((Ppos % 2) ? x * pow<Ppos - 1>(x) : pow<Ppos / 2>(x) * pow<Ppos / 2>(x))
             : 1 / pow<Ppos>(x);
}

template <>
PORTABLE_FORCEINLINE_FUNCTION constexpr Real pow<0>(const Real x) {
  return 1.0;
}

template <>
PORTABLE_FORCEINLINE_FUNCTION constexpr Real pow<1>(const Real x) {
  return x;
}

template <int P>
PORTABLE_FORCEINLINE_FUNCTION constexpr auto ipow10() {
  return pow<P>(10);
}

PORTABLE_FORCEINLINE_FUNCTION auto pow10(const Real x) {
  constexpr Real ln10 = 2.30258509299405e+00;
  return std::exp(ln10 * x);
}

} // namespace math_utils
} // namespace singularity

#endif // SINGULARITY_EOS_BASE_MATH_UTILS_HPP_

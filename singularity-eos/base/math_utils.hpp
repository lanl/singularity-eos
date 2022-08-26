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
 * Generic is recursive, but we specialize a few of them strategically
 * A generalization would be to treat even powers specially
 */
template <size_t P>
PORTABLE_FORCEINLINE_FUNCTION constexpr auto pow(const Real x) {
  return x * pow<P - 1>(x);
}

template <>
PORTABLE_FORCEINLINE_FUNCTION constexpr auto pow<0>(const Real x) {
  return 1.0;
}

template <>
PORTABLE_FORCEINLINE_FUNCTION constexpr auto pow<1>(const Real x) {
  return x;
}

template <>
PORTABLE_FORCEINLINE_FUNCTION constexpr auto pow<2>(const Real x) {
  return x * x;
}

template <>
PORTABLE_FORCEINLINE_FUNCTION constexpr auto pow<4>(const Real x) {
  return pow<2>(x) * pow<2>(x);
}

} // namespace math_utils
} // namespace singularity

#endif // SINGULARITY_EOS_BASE_MATH_UTILS_HPP_

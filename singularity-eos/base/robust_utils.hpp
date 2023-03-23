//------------------------------------------------------------------------------
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
//------------------------------------------------------------------------------

#ifndef SINGULARITY_EOS_BASE_ROBUST_UTILS_HPP_
#define SINGULARITY_EOS_BASE_ROBUST_UTILS_HPP_

#include <limits>
#include <ports-of-call/portability.hpp>

namespace singularity {
namespace robust {

template <typename T = Real>
PORTABLE_FORCEINLINE_FUNCTION constexpr auto SMALL() {
  return 10 * std::numeric_limits<T>::min();
}

template <typename T = Real>
PORTABLE_FORCEINLINE_FUNCTION constexpr auto EPS() {
  return 10 * std::numeric_limits<T>::epsilon();
}

template <typename T>
PORTABLE_FORCEINLINE_FUNCTION auto make_positive(const T val) {
  return std::max(val, EPS<T>());
}

PORTABLE_FORCEINLINE_FUNCTION
Real make_bounded(const Real val, const Real vmin, const Real vmax) {
  return std::min(std::max(val, vmin + EPS()), vmax * (1.0 - EPS()));
}

template <typename T>
PORTABLE_FORCEINLINE_FUNCTION int sgn(const T &val) {
  return (T(0) <= val) - (val < T(0));
}

template <typename A, typename B>
PORTABLE_FORCEINLINE_FUNCTION auto ratio(const A &a, const B &b) {
  return a / (b + sgn(b) * SMALL<B>());
}

} // namespace robust
} // namespace singularity

#endif // SINGULARITY_EOS_BASE_ROBUST_UTILS_HPP_

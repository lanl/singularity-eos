//------------------------------------------------------------------------------
// © 2021-2025. Triad National Security, LLC. All rights reserved.  This
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

#ifndef SINGULARITY_EOS_BASE_ERROR_UTILS_HPP_
#define SINGULARITY_EOS_BASE_ERROR_UTILS_HPP_

#include <cmath>
#include <type_traits>

#include <ports-of-call/portability.hpp>

namespace singularity {
namespace error_utils {

using PortsOfCall::printf;

constexpr double _NORMAL_FACTOR = 1.0e10;

struct is_normal_or_zero {
  template <typename valT>
  constexpr bool PORTABLE_FORCEINLINE_FUNCTION operator()(valT value) const {
    static_assert(std::is_floating_point<valT>::value);
    return (value == valT{0}) ||
           (std::isnormal(_NORMAL_FACTOR * value) && std::isnormal(value));
  }
};

struct is_strictly_positive {
  template <typename valT>
  constexpr bool PORTABLE_FORCEINLINE_FUNCTION operator()(valT value) const {
    static_assert(std::is_arithmetic<valT>::value);
    return value > valT{0};
  }
};

struct is_non_negative {
  template <typename valT>
  constexpr bool PORTABLE_FORCEINLINE_FUNCTION operator()(valT value) const {
    static_assert(std::is_arithmetic<valT>::value);
    return value >= valT{0};
  }
};

// Checks whether a value obeys some sort of provided condition. If not, returns true and
// prints the provided error message, variable name, and value (but does not abort!)
template <typename valT, typename condT, typename nameT>
PORTABLE_FORCEINLINE_FUNCTION bool violates_condition(valT &&value, condT &&condition,
                                                      nameT &&var_name) {
  constexpr bool good = condition(std::forward<valT>(value));
  if (!good) {
    printf("### ERROR: Bad singularity-eos value\n  Var:   %s\n  Value: %.15e\n",
           var_name, value);
  }
  return !good;
}

// Create specialized functions using the above conditions
template <typename valT, typename nameT>
PORTABLE_FORCEINLINE_FUNCTION bool bad_value(valT &&value, nameT &&var_name) {
  return violates_condition(std::forward<valT>(value), is_normal_or_zero{},
                            std::forward<nameT>(var_name));
}
template <typename valT, typename nameT>
PORTABLE_FORCEINLINE_FUNCTION bool non_positive_value(valT &&value, nameT &&var_name) {
  return violates_condition(std::forward<valT>(value), is_strictly_positive{},
                            std::forward<nameT>(var_name));
}
template <typename valT, typename nameT>
PORTABLE_FORCEINLINE_FUNCTION bool negative_value(valT &&value, nameT &&var_name) {
  return violates_condition(std::forward<valT>(value), is_strictly_positive{},
                            std::forward<nameT>(var_name));
}

} // namespace error_utils
} // namespace singularity

#endif // #ifndef SINGULARITY_EOS_BASE_ERROR_UTILS_HPP_

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

#ifndef _SINGULARITY_EOS_EOS_SINGULARITY_EOS_INIT_UTILS_HPP_
#define _SINGULARITY_EOS_EOS_SINGULARITY_EOS_INIT_UTILS_HPP_

#include <type_traits>
#include <utility>

#include <ports-of-call/portability.hpp>
#include <singularity-eos/base/variadic_utils.hpp>
#include <singularity-eos/eos/eos.hpp>

namespace singularity {
namespace init_utils_impl {
// true if T is one of the alternatives of the variant V
template <typename T, typename V>
struct in_variant : std::false_type {};
template <typename T, typename... Ts>
struct in_variant<T, Variant<Ts...>> : variadic_utils::contains<T, Ts...> {};

// Convert to EOS if T is an alternative of the EOS variant, otherwise
// abort at runtime. This lets the fortran backend compile against a
// custom SINGULARITY_VARIANT that does not contain every model and
// modifier combination the init functions can produce.
template <typename T>
inline EOS toEOS(T &&eos) {
  if constexpr (in_variant<std::decay_t<T>, EOS>::value) {
    return EOS(std::forward<T>(eos));
  } else {
    PORTABLE_ALWAYS_THROW_OR_ABORT("Requested EOS model/modifier combination is not in "
                                   "the singularity::EOS variant");
    return EOS();
  }
}
} // namespace init_utils_impl

// TODO: Replace these with the new Modify method in EOSBuilder.
// NOTE: The new EOSBuilder machinery will likely be slower than these.
template <typename T>
inline EOS applyShiftAndScale(T &&eos, bool scaled, bool shifted, Real scale,
                              Real shift) {
  if (shifted && scaled) {
    ShiftedEOS<T> a(std::forward<T>(eos), shift);
    ScaledEOS<ShiftedEOS<T>> b(std::move(a), scale);
    return init_utils_impl::toEOS(std::move(b));
  }
  if (shifted) {
    return init_utils_impl::toEOS(ShiftedEOS<T>(std::forward<T>(eos), shift));
  }
  if (scaled) {
    return init_utils_impl::toEOS(ScaledEOS<T>(std::forward<T>(eos), scale));
  }
  return init_utils_impl::toEOS(std::forward<T>(eos));
}

template <typename T, template <typename> class W, typename... ARGS>
inline EOS applyWrappedShiftAndScale(T &&eos, bool scaled, bool shifted, Real scale,
                                     Real shift, ARGS... args) {
  if (shifted && scaled) {
    ShiftedEOS<T> a(std::forward<T>(eos), shift);
    ScaledEOS<ShiftedEOS<T>> b(std::move(a), scale);
    W<ScaledEOS<ShiftedEOS<T>>> c(std::move(b), args...);
    return init_utils_impl::toEOS(std::move(c));
  }
  if (shifted) {
    ShiftedEOS<T> sh_eos(std::forward<T>(eos), shift);
    return init_utils_impl::toEOS(W<ShiftedEOS<T>>(std::move(sh_eos), args...));
  }
  if (scaled) {
    ScaledEOS<T> sc_eos(std::forward<T>(eos), scale);
    return init_utils_impl::toEOS(W<ScaledEOS<T>>(std::move(sc_eos), args...));
  }
  return init_utils_impl::toEOS(W<T>(std::forward<T>(eos), args...));
}

template <typename T>
inline EOS applyShiftAndScaleAndBilinearRamp(T &&eos, bool scaled, bool shifted,
                                             bool ramped, Real scale, Real shift, Real r0,
                                             Real a, Real b, Real c) {
  if (ramped) {
    return applyWrappedShiftAndScale<T, BilinearRampEOS>(
        std::forward<T>(eos), scaled, shifted, scale, shift, r0, a, b, c);
  } else {
    return applyShiftAndScale(std::forward<T>(eos), scaled, shifted, scale, shift);
  }
}

// apply everything but ramp in order to possibly calculate the
// SAP ramp parameters from p-alhpa ramp parameters
#define SGAPPLYMODSIMPLE(A)                                                              \
  applyShiftAndScale(A, enabled[0] == 1, enabled[1] == 1, vals[0], vals[1])

#define SGAPPLYMOD(A)                                                                    \
  applyShiftAndScaleAndBilinearRamp(A, enabled[0] == 1, enabled[1] == 1,                 \
                                    enabled[2] == 1 || enabled[3] == 1, vals[0],         \
                                    vals[1], vals[2], vals[3], vals[4], vals[5])

extern int def_en[4];
extern double def_v[6];

} // namespace singularity
#endif // _SINGULARITY_EOS_EOS_SINGULARITY_EOS_INIT_UTILS_HPP_

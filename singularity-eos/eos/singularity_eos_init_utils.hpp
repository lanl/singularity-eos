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

#include <ports-of-call/portability.hpp>
#include <singularity-eos/eos/eos.hpp>

namespace singularity {
// TODO: Replace these with the new Modify method in EOSBuilder.
// NOTE: The new EOSBuilder machinery will likely be slower than these.
template <typename T>
inline EOS applyShiftAndScale(T &&eos, bool scaled, bool shifted, Real scale,
                              Real shift) {
  if (shifted && scaled) {
    ShiftedEOS<T> a(std::forward<T>(eos), shift);
    ScaledEOS<ShiftedEOS<T>> b(std::move(a), scale);
    return b;
  }
  if (shifted) {
    return ShiftedEOS<T>(std::forward<T>(eos), shift);
  }
  if (scaled) {
    return ScaledEOS<T>(std::forward<T>(eos), scale);
  }
  return eos;
}

template <typename T, template <typename> class W, typename... ARGS>
inline EOS applyWrappedShiftAndScale(T &&eos, bool scaled, bool shifted, Real scale,
                                     Real shift, ARGS... args) {
  if (shifted && scaled) {
    ShiftedEOS<T> a(std::forward<T>(eos), shift);
    ScaledEOS<ShiftedEOS<T>> b(std::move(a), scale);
    W<ScaledEOS<ShiftedEOS<T>>> c(std::move(b), args...);
    return c;
  }
  if (shifted) {
    ShiftedEOS<T> sh_eos(std::forward<T>(eos), shift);
    return W<ShiftedEOS<T>>(std::move(sh_eos), args...);
  }
  if (scaled) {
    ScaledEOS<T> sc_eos(std::forward<T>(eos), scale);
    return W<ScaledEOS<T>>(std::move(sc_eos), args...);
  }
  return W<T>(std::forward<T>(eos), args...);
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

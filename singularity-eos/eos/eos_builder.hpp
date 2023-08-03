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

#ifndef _SINGULARITY_EOS_EOS_EOS_BUILDER_HPP_
#define _SINGULARITY_EOS_EOS_EOS_BUILDER_HPP_

#include <map>
#include <mpark/variant.hpp>
#include <ports-of-call/portability.hpp>
#include <singularity-eos/eos/eos.hpp>
#include <string>
#include <unordered_set>

namespace singularity {

namespace EOSBuilder {

// TODO(JMM): This could be strings? Or at least use a macro so we
// only have to write these down once?
// strings might allow us to create these structrs more easily
// in the host code.
enum class EOSType {
  IdealGas,
  Gruneisen,
  JWL,
  DavisReactants,
  DavisProducts
#ifdef SINGULARITY_USE_SPINER_WITH_HDF5
  ,
  SpinerEOSDependsRhoT,
  SpinerEOSDependsRhoSie,
  StellarCollapse
#endif
};
enum class EOSModifier { Scaled, Shifted, Relativistic, UnitSystem, BilinearRamp };

// evil type erasure
using param_t = mpark::variant<bool, int, Real, std::string>;
using params_t = std::map<std::string, param_t>;
using modifiers_t = std::map<EOSModifier, params_t>;
const params_t NO_PARAMS;

// TODO: Extend if needed
const std::unordered_set<EOSType> modifiable({EOSType::IdealGas
#ifdef SINGULARITY_USE_SPINER_WITH_HDF5
                                              ,
                                              EOSType::SpinerEOSDependsRhoT,
                                              EOSType::SpinerEOSDependsRhoSie,
                                              EOSType::StellarCollapse
#endif
});
bool isModifiable(EOSType t);
EOS buildEOS(EOSType type, params_t base_params, modifiers_t modifiers);
inline auto buildEOS(EOSType type, params_t base_params) {
  modifiers_t modifiers;
  return buildEOS(type, base_params, modifiers);
}

template <typename T>
EOS makeUnitSystem(T &&eos, bool use_length_time, Real rho_unit, Real sie_unit,
                   Real temp_unit, Real time_unit, Real mass_unit, Real length_unit) {
  if (use_length_time) {
    return UnitSystem<T>(std::move(eos), eos_units_init::length_time_units_init_tag,
                         time_unit, mass_unit, length_unit, temp_unit);
  } else {
    return UnitSystem<T>(std::move(eos), eos_units_init::thermal_units_init_tag, rho_unit,
                         sie_unit, temp_unit);
  }
}

template <typename T>
EOS makeRelativistic(T &&eos, Real cl) {
  return RelativisticEOS<T>(std::move(eos), cl);
}

template <typename T>
EOS makeBilinearRamp(T &&eos, Real r0, Real a, Real b, Real c) {
  return BilinearRampEOS<T>(std::forward<T>(eos), r0, a, b, c);
}

template <typename T>
EOS applyShiftAndScale(T &&eos, bool scaled, bool shifted, Real scale, Real shift) {
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
EOS applyWrappedShiftAndScale(T &&eos, bool scaled, bool shifted, Real scale, Real shift,
                              ARGS... args) {
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
EOS applyShiftAndScaleAndBilinearRamp(T &&eos, bool scaled, bool shifted, bool ramped,
                                      Real scale, Real shift, Real r0, Real a, Real b,
                                      Real c) {
  if (ramped) {
    return applyWrappedShiftAndScale<T, BilinearRampEOS>(
        std::forward<T>(eos), scaled, shifted, scale, shift, r0, a, b, c);
  } else {
    return applyShiftAndScale(std::forward<T>(eos), scaled, shifted, scale, shift);
  }
}
} // namespace EOSBuilder

} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_EOS_BUILDER_HPP_

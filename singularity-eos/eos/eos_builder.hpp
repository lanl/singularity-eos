//------------------------------------------------------------------------------
// © 2021. Triad National Security, LLC. All rights reserved.  This
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
#include <string>
#include <unordered_set>
#include <mpark/variant.hpp>
#include <ports-of-call/portability.hpp>
#include <singularity-eos/eos/eos.hpp>

namespace singularity {

namespace EOSBuilder {

  // TODO(JMM): This could be strings? Or at least use a macro so we
  // only have to write these down once?
  // strings might allow us to create these structrs more easily
  // in the host code.
  enum class EOSType {
    IdealGas
    , Gruneisen
    , JWL
    , DavisReactants
    , DavisProducts
#ifdef SPINER_USE_HDF
    , SpinerEOSDependsRhoT
    , SpinerEOSDependsRhoSie
    , StellarCollapse
#endif
  };
  enum class EOSModifier { Scaled, Shifted, Relativistic };

  // evil type erasure
  using param_t = mpark::variant<bool,int,Real,std::string>;
  using params_t = std::map<std::string,param_t>;
  using modifiers_t = std::map<EOSModifier,params_t>;
  const params_t NO_PARAMS;

  const std::unordered_set<EOSType>
  modifiable({EOSType::IdealGas
#ifdef SPINER_USE_HDF
	, EOSType::SpinerEOSDependsRhoT
	, EOSType::SpinerEOSDependsRhoSie
  , EOSType::StellarCollapse
#endif
	});
  bool isModifiable(EOSType t);
  EOS buildEOS(EOSType type,
				 params_t base_params,
				 modifiers_t modifiers);
  inline auto buildEOS(EOSType type, params_t base_params) {
    modifiers_t modifiers;
    return buildEOS(type,base_params,modifiers);
  }

  template<typename T>
  EOS makeRelativistic(T&& eos, Real cl) {
    return RelativisticEOS<T>(std::move(eos), cl);
  }

  template<typename T>
  EOS applyScaleAndShift(T&& eos,
			 bool scaled, bool shifted,
			 Real scale, Real shift) {
    if (shifted && scaled) {
      ScaledEOS<T> a(std::move(eos),scale);
      ShiftedEOS<ScaledEOS<T>> b(std::move(a),shift);
      return b;
    }
    if (shifted) {
      return ShiftedEOS<T>(std::move(eos),shift);
    }
    if (scaled) {
      return ScaledEOS<T>(std::move(eos),scale);
    }
    return eos;
  }

  template<typename T>
  EOS applyShiftAndScale(T&& eos,
                         bool scaled, bool shifted,
                         Real scale, Real shift) {
    if (shifted && scaled) {
      ShiftedEOS<T> a(std::forward<T>(eos),shift);
      ScaledEOS<ShiftedEOS<T>> b(std::move(a),scale);
      return b;
    }
    if (shifted) {
      return ShiftedEOS<T>(std::forward<T>(eos),shift);
    }
    if (scaled) {
      return ScaledEOS<T>(std::forward<T>(eos),scale);
    }
    return eos;
  }
} // namespace EOSBuilder

} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_EOS_BUILDER_HPP_

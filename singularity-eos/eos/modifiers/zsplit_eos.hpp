//------------------------------------------------------------------------------
// © 2024. Triad National Security, LLC. All rights reserved.  This
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

#ifndef _SINGULARITY_EOS_EOS_ZSPLIT_EOS_
#define _SINGULARITY_EOS_EOS_ZSPLIT_EOS_

#include <cmath>
#include <utility>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_errors.hpp>
#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/eos_error.hpp>
#include <singularity-eos/base/robust_utils.hpp>
#include <singularity-eos/eos/eos_base.hpp>

namespace singularity {

using namespace eos_base;

enum class ZSplitComponent { Electrons, Ions };

// TODO(JMM): Separate type? Or runtime?
template <ZSplitComponent ztype, typename T>
class ZSplit : public EosBase<ZSplit<ztype, T>> {
 public:
  using BaseType = T;
  static constexpr const ZSplitComponent zsplit_type = ztype;
  SG_ADD_BASE_CLASS_USINGS(ZSplit<ztype, T>);

  static std::string EosType() {
    std::string ts = (ztype == ZSplitComponent::Electrons) ? "Electrons, " : "Ions, ";
    return std::string("ZSplit<") + ts + T::EosType() + std::string(">");
  }
  static std::string EosPyType() {
    std::string ts = (ztype == ZSplitComponent::Electrons) ? "Electrons" : "Ions";
    return std::string("ZSplit") + ts + T::EosPyType();
  }

  ZSplit() = default;
  PORTABLE_INLINE_FUNCTION
  ZSplit(T &&t) : t_(std::forward<T>(t)) {}

  PORTABLE_INLINE_FUNCTION void CheckParams() const { t_.CheckParams(); }
  auto GetOnDevice() { return ZSplit<ztype, T>(t_.GetOnDevice()); }
  inline void Finalize() { t_.Finalize(); }

  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    const Real iscale = GetInvScale_(lambda);
    return t_.TemperatureFromDensityInternalEnergy(rho, sie * iscale, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temperature, Indexer_t &&lambda = nullptr) const {
    const Real scale = GetScale_(lambda);
    const Real et = t_.InternalEnergyFromDensityTemperature(rho, temperature, lambda);
    return scale * et;
  }

  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real PressureFromDensityTemperature(
      const Real rho, const Real temperature, Indexer_t &&lambda = nullptr) const {
    const Real scale = GetScale_(lambda);
    const Real Pt = t_.PressureFromDensityTemperature(rho, temperature, lambda);
    return scale * Pt;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real iscale = GetInvScale_(lambda);
    const Real scale = GetScale_(lambda);
    const Real Pt = t_.PressureFromDensityInternalEnergy(rho, iscale * sie, lambda);
    return scale * Pt;
  }

  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real
  MinInternalEnergyFromDensity(const Real rho, Indexer_t &&lambda = nullptr) const {
    const Real scale = GetScale_(lambda);
    const Real et = t_.MinInternalEnergyFromDensity(rho, lambda);
    return scale * et;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real
  EntropyFromDensityTemperature(const Real rho, const Real temperature,
                                Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real scale = GetScale_(lambda);
    return scale * t_.EntropyFromDensityTemperature(rho, temperature, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real EntropyFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    const Real scale = GetScale_(lambda);
    const Real iscale = GetInvScale_(lambda);
    return scale * t_.EntropyFromDensityInternalEnergy(rho, iscale * sie, lambda);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real scale = GetScale_(lambda);
    return scale * t_.SpecificHeatFromDensityTemperature(rho, temperature, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real scale = GetScale_(lambda);
    const Real iscale = GetInvScale_(lambda);
    return scale * t_.SpecificHeatFromDensityInternalEnergy(rho, iscale * sie, lambda);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real scale = GetScale_(lambda);
    return scale * t_.BulkModulusFromDensityTemperature(rho, temperature, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real scale = GetScale_(lambda);
    const Real iscale = GetInvScale_(lambda);
    return scale * t_.BulkModulusFromDensityInternalEnergy(rho, iscale * sie, lambda);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    // (1/rho) (dP/de). Scale cancels in P and e
    return t_.GruneisenParameterFromDensityTemperature(rho, temperature, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real iscale = GetInvScale_(lambda);
    return t_.GruneisenParameterFromDensityInternalEnergy(rho, iscale * sie, lambda);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  FillEos(Real &rho, Real &temp, Real &energy, Real &press, Real &cv, Real &bmod,
          const unsigned long output,
          Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real scale = GetScale_(lambda);
    const Real iscale = GetInvScale_(lambda);

    const bool sie_output = (output & thermalqs::specific_internal_energy);
    const bool sie_input = !sie_output;

    const bool press_output = (output & thermalqs::pressure);
    const bool press_input = !press_output;

    const bool cv_output = (output & thermalqs::specific_heat);
    const bool bmod_output = (output & thermalqs::bulk_modulus);

    if (sie_input) energy *= iscale;
    if (press_input) press *= iscale;

    t_.FillEos(rho, temp, energy, press, cv, output, lambda);

    // do it to undo if it was input. do it because you need to if
    // output
    energy *= scale;
    press *= scale;
    // These aren't input
    if (cv_output) cv *= scale;
    if (bmod_output) bmod *= scale;
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                         Real &bmod, Real &dpde, Real &dvdt,
                         Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real scale = GetScale_(lambda);
    const Real iscale = GetInvScale_(lambda);
    t_.ValuesAtReferenceState(rho, temp, sie, press, cv, bmod, dpde, dvdt);
    sie *= scale;
    press *= scale;
    cv *= scale;
    bmod *= scale;
    // dvdt, dpde unchanged
  }

  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return NL + t_.nlambda(); }
  static constexpr unsigned long PreferredInput() { return T::PreferredInput(); }
  static inline unsigned long scratch_size(std::string method, unsigned int nelements) {
    return T::scratch_size(method, nelements);
  }
  static inline unsigned long max_scratch_size(unsigned int nelements) {
    return T::max_scratch_size(nelements);
  }
  PORTABLE_FUNCTION void PrintParams() const { t_.PrintParams(); }

  PORTABLE_FORCEINLINE_FUNCTION Real MinimumDensity() const {
    return t_.MinimumDensity();
  }
  PORTABLE_FORCEINLINE_FUNCTION Real MinimumTemperature() const {
    return t_.MinimumTemperature();
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real MeanAtomicMassFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return t_.MeanAtomicMassFromDensityTemperature(rho, temperature, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real MeanAtomicNumberFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return t_.MeanAtomicNumberFromDensityTemperature(rho, temperature, lambda);
  }

  SG_ADD_MODIFIER_METHODS(T, t_);
  SG_ADD_MODIFIER_MEAN_METHODS(t_)

 private:
  static constexpr const int NL = 1;
  template <typename Indexer_t>
  PORTABLE_FORCEINLINE_FUNCTION Real GetIonizationState_(const Indexer_t &&lambda) const {
    return std::max(0.0, lambda[t_.nlambda()]);
  }
  // TODO(JMM): Runtime?
  template <typename Indexer_t>
  PORTABLE_FORCEINLINE_FUNCTION Real GetScale_(const Indexer_t &&lambda) const {
    Real Z = GetIonizationState_(lambda);
    if constexpr (ztype == ZSplitComponent::Electrons) {
      return robust::ratio(Z, Z + 1);
    } else {
      return robust::ratio(1.0, Z + 1);
    }
  }
  template <typename Indexer_t>
  PORTABLE_FORCEINLINE_FUNCTION Real GetInvScale_(const Indexer_t &&lambda) const {
    Real Z = GetIonizationState_(lambda);
    if constexpr (ztype == ZSplitComponent::Electrons) {
      return robust::ratio(Z + 1, Z);
    } else {
      return Z + 1;
    }
  }
  T t_;
};

template <typename T>
using ZSplitE = ZSplit<ZSplitComponent::Electrons, T>;

template <typename T>
using ZSplitI = ZSplit<ZSplitComponent::Ions, T>;

} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_ZSPLIT_EOS_

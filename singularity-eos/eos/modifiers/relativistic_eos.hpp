//------------------------------------------------------------------------------
// © 2021-2024. Triad National Security, LLC. All rights reserved.  This
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

#ifndef _SINGULARITY_EOS_EOS_RELATIVISTIC_EOS_
#define _SINGULARITY_EOS_EOS_RELATIVISTIC_EOS_

#include "stdio.h"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <utility>

#include <ports-of-call/portability.hpp>
#include <singularity-eos/base/eos_error.hpp>
#include <singularity-eos/base/robust_utils.hpp>
#include <singularity-eos/eos/eos_base.hpp>

namespace singularity {

using namespace eos_base;

template <typename T>
class RelativisticEOS : public EosBase<RelativisticEOS<T>> {
 public:
  // Generic functions provided by the base class. These contain
  // e.g. the vector overloads that use the scalar versions declared
  // here We explicitly list, rather than using the macro because we
  // overload some methods.

  // TODO(JMM): The modifier EOS's should probably call the specific
  // sub-functions of the class they modify so that they can leverage,
  // e.g., an especially performant or special version of these
  using EosBase<RelativisticEOS<T>>::TemperatureFromDensityInternalEnergy;
  using EosBase<RelativisticEOS<T>>::InternalEnergyFromDensityTemperature;
  using EosBase<RelativisticEOS<T>>::PressureFromDensityTemperature;
  using EosBase<RelativisticEOS<T>>::PressureFromDensityInternalEnergy;
  using EosBase<RelativisticEOS<T>>::MinInternalEnergyFromDensity;
  using EosBase<RelativisticEOS<T>>::EntropyFromDensityTemperature;
  using EosBase<RelativisticEOS<T>>::EntropyFromDensityInternalEnergy;
  using EosBase<RelativisticEOS<T>>::SpecificHeatFromDensityTemperature;
  using EosBase<RelativisticEOS<T>>::SpecificHeatFromDensityInternalEnergy;
  using EosBase<RelativisticEOS<T>>::BulkModulusFromDensityTemperature;
  using EosBase<RelativisticEOS<T>>::BulkModulusFromDensityInternalEnergy;
  using EosBase<RelativisticEOS<T>>::GruneisenParamFromDensityTemperature;
  using EosBase<RelativisticEOS<T>>::GruneisenParamFromDensityInternalEnergy;
  using EosBase<RelativisticEOS<T>>::FillEos;
  using EosBase<RelativisticEOS<T>>::SerializedSizeInBytes;
  using EosBase<RelativisticEOS<T>>::Serialize;
  using EosBase<RelativisticEOS<T>>::DeSerialize;

  using BaseType = T;

  // give me std::format or fmt::format...
  static std::string EosType() {
    return std::string("RelativisticEOS<") + T::EosType() + std::string(">");
  }

  static std::string EosPyType() { return std::string("Relativistic") + T::EosPyType(); }

  // move semantics ensures dynamic memory comes along for the ride
  RelativisticEOS(T &&t, const Real cl)
      : t_(std::forward<T>(t)), cl_(cl) // speed of light, units arbitrary
        ,
        cl2_(cl * cl) // speed of light squared
  {}
  RelativisticEOS() = default;

  auto GetOnDevice() { return RelativisticEOS<T>(t_.GetOnDevice(), cl_); }
  inline void Finalize() { t_.Finalize(); }

  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    return t_.TemperatureFromDensityInternalEnergy(rho, sie, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temperature, Indexer_t &&lambda = nullptr) const {
    return t_.InternalEnergyFromDensityTemperature(rho, temperature, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real PressureFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    return t_.PressureFromDensityInternalEnergy(rho, sie, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real
  MinInternalEnergyFromDensity(const Real rho, Indexer_t &&lambda = nullptr) const {
    return t_.MinInternalEnergyFromDensity(rho, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real EntropyFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    return t_.EntropyFromDensityInternalEnergy(rho, sie, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    return t_.SpecificHeatFromDensityInternalEnergy(rho, sie, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    Real P = PressureFromDensityInternalEnergy(rho, sie, lambda);
    Real h = cl2_ + sie + robust::ratio(P, rho);
    Real bmod = t_.BulkModulusFromDensityInternalEnergy(rho, sie, lambda);
    return std::max(0.0, robust::ratio(bmod, std::abs(h)));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    return t_.GruneisenParamFromDensityInternalEnergy(rho, sie, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real PressureFromDensityTemperature(
      const Real rho, const Real temperature, Indexer_t &&lambda = nullptr) const {
    return t_.PressureFromDensityTemperature(rho, temperature, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real EntropyFromDensityTemperature(
      const Real rho, const Real temperature, Indexer_t &&lambda = nullptr) const {
    return t_.EntropyFromDensityTemperature(rho, temperature, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real SpecificHeatFromDensityTemperature(
      const Real rho, const Real temperature, Indexer_t &&lambda = nullptr) const {
    return t_.SpecificHeatFromDensityTemperature(rho, temperature, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real BulkModulusFromDensityTemperature(
      const Real rho, const Real temperature, Indexer_t &&lambda = nullptr) const {
    Real P = PressureFromDensityTemperature(rho, temperature, lambda);
    Real sie = InternalEnergyFromDensityTemperature(rho, temperature, lambda);
    Real h = cl2_ + sie + robust::ratio(P, rho);
    Real bmod = t_.BulkModulusFromDensityTemperature(rho, temperature, lambda);
    return std::max(0.0, robust::ratio(bmod, std::abs(h)));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temperature, Indexer_t &&lambda = nullptr) const {
    return t_.GruneisenParamFromDensityTemperature(rho, temperature, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION void FillEos(Real &rho, Real &temp, Real &energy, Real &press,
                                 Real &cv, Real &bmod, const unsigned long output,
                                 Indexer_t &&lambda = nullptr) const {
    t_.FillEos(rho, temp, energy, press, cv, bmod, output, lambda);
  }

  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return t_.nlambda(); }

  PORTABLE_FORCEINLINE_FUNCTION Real MinimumDensity() const {
    return t_.MinimumDensity();
  }
  PORTABLE_FORCEINLINE_FUNCTION Real MinimumTemperature() const {
    return t_.MinimumTemperature();
  }

  static constexpr unsigned long PreferredInput() { return T::PreferredInput(); }

  static inline unsigned long scratch_size(std::string method, unsigned int nelements) {
    return T::scratch_size(method, nelements);
  }

  static inline unsigned long max_scratch_size(unsigned int nelements) {
    return T::max_scratch_size(nelements);
  }

  PORTABLE_FUNCTION void PrintParams() const { t_.PrintParams(); }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION void
  DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                       Indexer_t &&lambda, Real &rho, Real &sie) const {
    t_.DensityEnergyFromPressureTemperature(press, temp, lambda, rho, sie);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION void ValuesAtReferenceState(Real &rho, Real &temp, Real &sie,
                                                Real &press, Real &cv, Real &bmod,
                                                Real &dpde, Real &dvdt,
                                                Indexer_t &&lambda = nullptr) const {
    t_.ValuesAtReferenceState(rho, temp, sie, press, cv, bmod, dpde, dvdt, lambda);
  }

  inline constexpr bool IsModified() const { return true; }

  inline constexpr T UnmodifyOnce() { return t_; }

  inline constexpr decltype(auto) GetUnmodifiedObject() {
    return t_.GetUnmodifiedObject();
  }

  std::size_t DynamicMemorySizeInBytes() const { return t_.DynamicMemorySizeInBytes(); }
  std::size_t DumpDynamicMemory(char *dst) const { return t_.DumpDynamicMemory(dst); }
  std::size_t SetDynamicMemory(char *src) { return t_.SetDynamicMemory(src); }

 private:
  T t_;
  Real cl_, cl2_;
};

} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_RELATIVISTIC_EOS_

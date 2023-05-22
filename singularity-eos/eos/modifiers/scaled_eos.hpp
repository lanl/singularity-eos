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

#ifndef _SINGULARITY_EOS_EOS_SCALED_EOS_
#define _SINGULARITY_EOS_EOS_SCALED_EOS_

#include "stdio.h"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <utility>

#include <ports-of-call/portability.hpp>
#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/eos_error.hpp>
#include <singularity-eos/eos/eos_base.hpp>

namespace singularity {

using namespace eos_base;

template <typename T>
class ScaledEOS : public EosBase<ScaledEOS<T>> {
 public:
  // Generic functions provided by the base class. These contain
  // e.g. the vector overloads that use the scalar versions declared
  // here We explicitly list, rather than using the macro because we
  // overload some methods.

  // TODO(JMM): The modifier EOS's should probably call the specific
  // sub-functions of the class they modify so that they can leverage,
  // e.g., an especially performant or special version of these
  using EosBase<ScaledEOS<T>>::TemperatureFromDensityInternalEnergy;
  using EosBase<ScaledEOS<T>>::InternalEnergyFromDensityTemperature;
  using EosBase<ScaledEOS<T>>::PressureFromDensityTemperature;
  using EosBase<ScaledEOS<T>>::PressureFromDensityInternalEnergy;
  using EosBase<ScaledEOS<T>>::EntropyFromDensityTemperature;
  using EosBase<ScaledEOS<T>>::EntropyFromDensityInternalEnergy;
  using EosBase<ScaledEOS<T>>::SpecificHeatFromDensityTemperature;
  using EosBase<ScaledEOS<T>>::SpecificHeatFromDensityInternalEnergy;
  using EosBase<ScaledEOS<T>>::BulkModulusFromDensityTemperature;
  using EosBase<ScaledEOS<T>>::BulkModulusFromDensityInternalEnergy;
  using EosBase<ScaledEOS<T>>::GruneisenParamFromDensityTemperature;
  using EosBase<ScaledEOS<T>>::GruneisenParamFromDensityInternalEnergy;
  using EosBase<ScaledEOS<T>>::PTofRE;
  using EosBase<ScaledEOS<T>>::FillEos;

  using BaseType = T;

  // give me std::format or fmt::format...
  static std::string EosType() {
    return std::string("ScaledEOS<") + T::EosType() + std::string(">");
  }

  static std::string EosPyType() { return std::string("Scaled") + T::EosPyType(); }

  // move semantics ensures dynamic memory comes along for the ride
  ScaledEOS(T &&t, const Real scale)
      : t_(std::forward<T>(t)), scale_(scale), inv_scale_(1. / scale) {}
  ScaledEOS() = default;

  auto GetOnDevice() { return ScaledEOS<T>(t_.GetOnDevice(), scale_); }
  inline void Finalize() { t_.Finalize(); }

  PORTABLE_FUNCTION
  Real TemperatureFromDensityInternalEnergy(const Real rho, const Real sie,
                                            Real *lambda = nullptr) const {
    return t_.TemperatureFromDensityInternalEnergy(scale_ * rho, inv_scale_ * sie,
                                                   lambda);
  }
  PORTABLE_FUNCTION
  Real InternalEnergyFromDensityTemperature(const Real rho, const Real temperature,
                                            Real *lambda = nullptr) const {
    Real energy = t_.InternalEnergyFromDensityTemperature(rho, temperature, lambda);
    return scale_ * energy;
  }
  PORTABLE_FUNCTION
  Real PressureFromDensityInternalEnergy(const Real rho, const Real sie,
                                         Real *lambda = nullptr) const {
    return t_.PressureFromDensityInternalEnergy(scale_ * rho, inv_scale_ * sie, lambda);
  }
  PORTABLE_FUNCTION
  Real EntropyFromDensityInternalEnergy(const Real rho, const Real sie,
                                        Real *lambda = nullptr) const {
    return scale_ *
           t_.EntropyFromDensityInternalEnergy(scale_ * rho, inv_scale_ * sie, lambda);
  }
  PORTABLE_FUNCTION
  Real SpecificHeatFromDensityInternalEnergy(const Real rho, const Real sie,
                                             Real *lambda = nullptr) const {
    return t_.SpecificHeatFromDensityInternalEnergy(scale_ * rho, inv_scale_ * sie,
                                                    lambda);
  }
  PORTABLE_FUNCTION
  Real BulkModulusFromDensityInternalEnergy(const Real rho, const Real sie,
                                            Real *lambda = nullptr) const {
    return t_.BulkModulusFromDensityInternalEnergy(scale_ * rho, inv_scale_ * sie,
                                                   lambda);
  }
  PORTABLE_FUNCTION
  Real GruneisenParamFromDensityInternalEnergy(const Real rho, const Real sie,
                                               Real *lambda = nullptr) const {
    return t_.GruneisenParamFromDensityInternalEnergy(scale_ * rho, inv_scale_ * sie,
                                                      lambda);
  }
  PORTABLE_FUNCTION
  Real PressureFromDensityTemperature(const Real rho, const Real temperature,
                                      Real *lambda = nullptr) const {
    return t_.PressureFromDensityTemperature(scale_ * rho, temperature, lambda);
  }
  PORTABLE_FUNCTION
  Real EntropyFromDensityTemperature(const Real rho, const Real temperature,
                                     Real *lambda = nullptr) const {
    return scale_ * t_.EntropyFromDensityTemperature(scale_ * rho, temperature, lambda);
  }
  PORTABLE_FUNCTION
  Real SpecificHeatFromDensityTemperature(const Real rho, const Real temperature,
                                          Real *lambda = nullptr) const {
    return t_.SpecificHeatFromDensityTemperature(scale_ * rho, temperature, lambda);
  }
  PORTABLE_FUNCTION
  Real BulkModulusFromDensityTemperature(const Real rho, const Real temperature,
                                         Real *lambda = nullptr) const {
    return t_.BulkModulusFromDensityTemperature(scale_ * rho, temperature, lambda);
  }
  PORTABLE_FUNCTION
  Real GruneisenParamFromDensityTemperature(const Real rho, const Real temperature,
                                            Real *lambda = nullptr) const {
    return t_.GruneisenParamFromDensityTemperature(scale_ * rho, temperature, lambda);
  }
  PORTABLE_FUNCTION
  void FillEos(Real &rho, Real &temp, Real &energy, Real &press, Real &cv, Real &bmod,
               const unsigned long output, Real *lambda = nullptr) const {
    Real srho, senergy;
    switch (t_.PreferredInput()) {
    case thermalqs::density | thermalqs::temperature:
      srho = scale_ * rho;
      t_.FillEos(srho, temp, energy, press, cv, bmod, output, lambda);
      energy = scale_ * energy;
      break;
    case thermalqs::density | thermalqs::specific_internal_energy:
      srho = scale_ * rho;
      senergy = inv_scale_ * energy;
      t_.FillEos(srho, temp, senergy, press, cv, bmod, output, lambda);
      break;
    default:
      EOS_ERROR("Didn't find a valid input for ScaledEOS::FillEOS\n");
    }
  }

  PORTABLE_FUNCTION
  void ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                              Real &bmod, Real &dpde, Real &dvdt,
                              Real *lambda = nullptr) const {
    t_.ValuesAtReferenceState(rho, temp, sie, press, cv, bmod, dpde, dvdt, lambda);
    rho *= inv_scale_;
    sie *= scale_;
  }

  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return t_.nlambda(); }

  static constexpr unsigned long PreferredInput() { return T::PreferredInput(); }

  static inline unsigned long scratch_size(std::string method, unsigned int nelements) {
    return T::scratch_size(method, nelements);
  }

  static inline unsigned long max_scratch_size(unsigned int nelements) {
    return T::max_scratch_size(nelements);
  }

  PORTABLE_FUNCTION void PrintParams() const {
    t_.PrintParams();
    printf("scaling_ratio = %f\n", scale_);
  }
  PORTABLE_FUNCTION
  void DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                            Real *lambda, Real &rho, Real &sie) const {
    t_.DensityEnergyFromPressureTemperature(press, temp, lambda, rho, sie);
    rho = rho * inv_scale_;
    sie = sie * scale_;
  }

  PORTABLE_FORCEINLINE_FUNCTION Real MinimumDensity() const {
    return inv_scale_ * t_.MinimumDensity();
  }
  PORTABLE_FORCEINLINE_FUNCTION Real MinimumTemperature() const {
    return t_.MinimumTemperature();
  }

  inline constexpr bool IsModified() const { return true; }

  inline constexpr T UnmodifyOnce() { return t_; }

  inline constexpr decltype(auto) GetUnmodifiedObject() { return t_.GetUnmodifiedObject(); }

 private:
  T t_;
  double scale_;
  double inv_scale_;
};

} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_SCALED_EOS_

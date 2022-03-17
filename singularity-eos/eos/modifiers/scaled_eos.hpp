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

  PORTABLE_FUNCTION
  unsigned long PreferredInput() const { return t_.PreferredInput(); }

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

  // Vector functions that overload the scalar versions declared here.
  SG_ADD_BASE_CLASS_USINGS(ScaledEOS<T>)

 private:
  T t_;
  double scale_;
  double inv_scale_;
};

} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_SCALED_EOS_

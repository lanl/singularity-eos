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

#ifndef _SINGULARITY_EOS_EOS_SHIFTED_EOS_
#define _SINGULARITY_EOS_EOS_SHIFTED_EOS_

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
class ShiftedEOS : public EosBase<ShiftedEOS<T>> {
 public:
  // move semantics ensures dynamic memory comes along for the ride
  ShiftedEOS(T &&t, const Real shift) : t_(std::forward<T>(t)), shift_(shift) {}
  ShiftedEOS() = default;

  auto GetOnDevice() { return ShiftedEOS<T>(t_.GetOnDevice(), shift_); }
  inline void Finalize() { t_.Finalize(); }

  PORTABLE_FUNCTION
  Real TemperatureFromDensityInternalEnergy(const Real rho, const Real sie,
                                            Real *lambda = nullptr) const {
    return t_.TemperatureFromDensityInternalEnergy(rho, sie - shift_, lambda);
  }
  PORTABLE_FUNCTION
  Real InternalEnergyFromDensityTemperature(const Real rho, const Real temperature,
                                            Real *lambda = nullptr) const {
    Real energy = t_.InternalEnergyFromDensityTemperature(rho, temperature, lambda);
    return energy + shift_;
  }
  PORTABLE_FUNCTION
  Real PressureFromDensityInternalEnergy(const Real rho, const Real sie,
                                         Real *lambda = nullptr) const {
    return t_.PressureFromDensityInternalEnergy(rho, sie - shift_, lambda);
  }
  PORTABLE_FUNCTION
  Real SpecificHeatFromDensityInternalEnergy(const Real rho, const Real sie,
                                             Real *lambda = nullptr) const {
    return t_.SpecificHeatFromDensityInternalEnergy(rho, sie - shift_, lambda);
  }
  PORTABLE_FUNCTION
  Real BulkModulusFromDensityInternalEnergy(const Real rho, const Real sie,
                                            Real *lambda = nullptr) const {
    return t_.BulkModulusFromDensityInternalEnergy(rho, sie - shift_, lambda);
  }
  PORTABLE_FUNCTION
  Real GruneisenParamFromDensityInternalEnergy(const Real rho, const Real sie,
                                               Real *lambda = nullptr) const {
    return t_.GruneisenParamFromDensityInternalEnergy(rho, sie - shift_, lambda);
  }
  PORTABLE_FUNCTION
  Real PressureFromDensityTemperature(const Real rho, const Real temperature,
                                      Real *lambda = nullptr) const {
    return t_.PressureFromDensityInternalEnergy(rho, temperature, lambda);
  }
  PORTABLE_FUNCTION
  Real SpecificHeatFromDensityTemperature(const Real rho, const Real temperature,
                                          Real *lambda = nullptr) const {
    return t_.SpecificHeatFromDensityTemperature(rho, temperature, lambda);
  }
  PORTABLE_FUNCTION
  Real BulkModulusFromDensityTemperature(const Real rho, const Real temperature,
                                         Real *lambda = nullptr) const {
    return t_.BulkModulusFromDensityTemperature(rho, temperature, lambda);
  }
  PORTABLE_FUNCTION
  Real GruneisenParamFromDensityTemperature(const Real rho, const Real temperature,
                                            Real *lambda = nullptr) const {
    return t_.GruneisenParamFromDensityTemperature(rho, temperature, lambda);
  }
  PORTABLE_FUNCTION
  void FillEos(Real &rho, Real &temp, Real &energy, Real &press, Real &cv, Real &bmod,
               const unsigned long output, Real *lambda = nullptr) const {
    Real senergy;
    switch (t_.PreferredInput()) {
    case thermalqs::density | thermalqs::temperature:
      t_.FillEos(rho, temp, energy, press, cv, bmod, output, lambda);
      energy = energy + shift_;
      break;
    case thermalqs::density | thermalqs::specific_internal_energy:
      senergy = energy - shift_;
      t_.FillEos(rho, temp, senergy, press, cv, bmod, output, lambda);
      break;
    default:
      EOS_ERROR("Didn't find a valid input for ShiftedEOS::FillEOS\n");
    }
  }

  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return t_.nlambda(); }

  PORTABLE_FUNCTION
  unsigned long PreferredInput() const { return t_.PreferredInput(); }

  PORTABLE_FUNCTION void PrintParams() const {
    t_.PrintParams();
    printf("scaling_ratio = %f\n", shift_);
  }
  PORTABLE_FUNCTION
  void DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                            Real *lambda, Real &rho, Real &sie) const {
    t_.DensityEnergyFromPressureTemperature(press, temp, lambda, rho, sie);
    sie = sie + shift_;
  }

  // PORTABLE_FUNCTION
  // void PTofRE(const Real rho, const Real sie, Real *lambda, Real &press, Real &temp,
  //             Real &dpdr, Real &dpde, Real &dtdr, Real &dtde) const {
  //   t_.PTofRE(rho, sie - shift_, lambda, press, temp, dpdr, dpde, dtdr, dtde);
  // }
  PORTABLE_FUNCTION
  void ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                              Real &bmod, Real &dpde, Real &dvdt,
                              Real *lambda = nullptr) const {
    t_.ValuesAtReferenceState(rho, temp, sie, press, cv, bmod, dpde, dvdt, lambda);
    sie += shift_;
  }

  // Vector functions that overload the scalar versions declared here.
  SG_ADD_BASE_CLASS_USINGS(ShiftedEOS<T>)

 private:
  T t_;
  double shift_;
};

} // namespace singularity

#endif _SINGULARITY_EOS_EOS_SHIFTED_EOS_

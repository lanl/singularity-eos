//------------------------------------------------------------------------------
// © 2021-2022. Triad National Security, LLC. All rights reserved.  This
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

#ifndef _SINGULARITY_EOS_EOS_RAMPS_EOS_
#define _SINGULARITY_EOS_EOS_RAMPS_EOS_

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
class SAPRampEOS : public EosBase<SAPRampEOS<T>> {
 public:
  // move semantics ensures dynamic memory comes along for the ride
  SAPRampEOS(T &&t, const Real r0, const Real a, const Real b, const Real c) 
  : t_(std::forward<T>(t)),
    r0_(r0),
    a_(a),
    b_(b),
    c_(c),
    rmid_(r0*(a-b*c)/(a-b))
    {}
  SAPRampEOS() = default;

  auto GetOnDevice() { return SAPRampEOS<T>(t_.GetOnDevice(), r0_, a_, b_, c_); }
  inline void Finalize() { t_.Finalize(); }

  PORTABLE_INLINE_FUNCTION
  Real get_ramp_pressure(Real rho) const {
    const Real p_ramp {rho < r0_ ? 0.0                  :
                       rho < rmid_ ? a_*(rho/r0_ - 1.0) :
                                       b_*(rho/r0_ - c_)};
    return p_ramp;
  }

  PORTABLE_FUNCTION
  Real TemperatureFromDensityInternalEnergy(const Real rho, const Real sie,
                                            Real *lambda = nullptr) const {
    return t_.TemperatureFromDensityInternalEnergy(rho, sie, lambda);
  }
  PORTABLE_FUNCTION
  Real InternalEnergyFromDensityTemperature(const Real rho, const Real temperature,
                                            Real *lambda = nullptr) const {
    return t_.InternalEnergyFromDensityTemperature(rho, temperature, lambda);
  }
  PORTABLE_FUNCTION
  Real PressureFromDensityInternalEnergy(const Real rho, const Real sie,
                                         Real *lambda = nullptr) const {
    // ramp pressure
    const Real p_ramp {get_ramp_pressure(rho)};
    // pressure from eos
    const Real p_eos {t_.PressureFromDensityInternalEnergy(rho, sie, lambda)};
    // return max(p_ramp, p_eos)
    return p_eos < p_ramp ? p_ramp : p_eos;
  }
  PORTABLE_FUNCTION
  Real SpecificHeatFromDensityInternalEnergy(const Real rho, const Real sie,
                                             Real *lambda = nullptr) const {
    return t_.SpecificHeatFromDensityInternalEnergy(rho, sie, lambda);
  }
  PORTABLE_FUNCTION
  Real BulkModulusFromDensityInternalEnergy(const Real rho, const Real sie,
                                            Real *lambda = nullptr) const {
    return t_.BulkModulusFromDensityInternalEnergy(rho, sie, lambda);
  }
  PORTABLE_FUNCTION
  Real GruneisenParamFromDensityInternalEnergy(const Real rho, const Real sie,
                                               Real *lambda = nullptr) const {
    return t_.GruneisenParamFromDensityInternalEnergy(rho, sie, lambda);
  }
  PORTABLE_FUNCTION
  Real PressureFromDensityTemperature(const Real rho, const Real temperature,
                                      Real *lambda = nullptr) const {
    const Real p_ramp {get_ramp_pressure(rho)};
    const Real p_eos {t_.PressureFromDensityTemperature(rho, temperature, lambda)};
    return p_eos < p_ramp ? p_ramp : p_eos;
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
    // density must be an input
    assert(!(output & thermalqs::density));
    // output passed into internal filleos can't include pressure
    const unsigned long ramp_out = output & ~thermalqs::pressure;
    // if pressure is output, calculate it first
    if (output & thermalqs::pressure) {
      // maybe switch case on preferred input and check for not output of the input
      // for now sie lookup is prioritized
      if (!(output & thermalqs::specific_internal_energy)) {
	press = this->PressureFromDensityInternalEnergy(rho, energy, lambda);
      }
      else if (!(output & thermalqs::temperature)) {
	press = this->PressureFromDensityTemperature(rho, temp, lambda);
      }
    }
    // call internal filleos
    t_.FillEos(rho, temp, energy, press, cv, bmod, ramp_out, lambda);
  }

  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return t_.nlambda(); }

  PORTABLE_FUNCTION
  unsigned long PreferredInput() const { return t_.PreferredInput(); }

  PORTABLE_FUNCTION void PrintParams() const {
    t_.PrintParams();
    printf("r0=%e\na=%e\nb=%e\nc=%e\n", r0_, a_, b_, c_);
  }
  PORTABLE_FUNCTION
  void DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                            Real *lambda, Real &rho, Real &sie) const {
    t_.DensityEnergyFromPressureTemperature(press, temp, lambda, rho, sie);
  }

  PORTABLE_FUNCTION
  void ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                              Real &bmod, Real &dpde, Real &dvdt,
                              Real *lambda = nullptr) const {
    t_.ValuesAtReferenceState(rho, temp, sie, press, cv, bmod, dpde, dvdt, lambda);
  }

  // Vector functions that overload the scalar versions declared here.
  SG_ADD_BASE_CLASS_USINGS(SAPRampEOS<T>)

 private:
  T t_;
  Real r0_;
  Real a_;
  Real b_;
  Real c_;
  Real rmid_;
};

} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_RAMPS_EOS_

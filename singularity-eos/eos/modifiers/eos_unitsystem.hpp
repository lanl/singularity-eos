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

#ifndef _SINGULARITY_EOS_EOS_EOS_UNITSYSTEM_HPP_
#define _SINGULARITY_EOS_EOS_EOS_UNITSYSTEM_HPP_

#include "stdio.h"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <utility>

#include <ports-of-call/portability.hpp>
#include <singularity-eos/eos/base/constants.hpp>
#include <singularity-eos/eos/base/eos_error.hpp>

namespace singularity {

// tag dispatch for constructors for UnitSystem
namespace eos_units_init {
static struct ThermalUnitsInit {
} thermal_units_init_tag;
static struct LengthTimeUnitsInit {
} length_time_units_init_tag;
} // namespace eos_units_init

template <typename T>
class UnitSystem {
 public:
  // move semantics ensures dynamic memory comes along for the ride
  // TODO(JMM): Entropy unit needed?
  UnitSystem(T &&t, eos_units_init::ThermalUnitsInit, const Real rho_unit,
             const Real sie_unit, const Real temp_unit)
      : t_(std::forward<T>(t)), rho_unit_(rho_unit), sie_unit_(sie_unit),
        temp_unit_(temp_unit), press_unit_(rho_unit * sie_unit),
        inv_rho_unit_(1 / rho_unit) // inverses computed to avoid division at runtime
        ,
        inv_sie_unit_(1 / sie_unit), inv_temp_unit_(1 / temp_unit),
        inv_press_unit_(1 / press_unit_),
        inv_dpde_unit_(sie_unit / press_unit_) // thermo derivatives computed consistently
        ,
        inv_dvdt_unit_(rho_unit * temp_unit) // TODO(JMM): Is this convention weird?
        ,
        inv_dpdr_unit_(rho_unit / press_unit_), inv_dtdr_unit_(rho_unit / temp_unit),
        inv_dtde_unit_(sie_unit / temp_unit) // obviously this is also Cv
        ,
        inv_cv_unit_(temp_unit / sie_unit),
        inv_bmod_unit_(1 / press_unit_) {}
  UnitSystem(T &&t, eos_units_init::LengthTimeUnitsInit, const Real time_unit,
             const Real mass_unit, const Real length_unit, const Real temp_unit)
      : UnitSystem(std::forward<T>(t), eos_units_init::thermal_units_init_tag,
                   mass_unit / (length_unit * length_unit * length_unit),
                   length_unit * length_unit / (time_unit * time_unit), temp_unit) {}
  UnitSystem(T &&t, const Real rho_unit, const Real sie_unit, const Real temp_unit)
      : UnitSystem(std::forward<T>(t), eos_units_init::thermal_units_init_tag, rho_unit,
                   sie_unit, temp_unit) {}
  UnitSystem() = default;

  auto GetOnDevice() {
    return UnitSystem<T>(t_.GetOnDevice(), eos_units_init::thermal_units_init_tag,
                         rho_unit_, sie_unit_, temp_unit_);
  }
  inline void Finalize() { t_.Finalize(); }

  PORTABLE_FUNCTION
  Real TemperatureFromDensityInternalEnergy(const Real rho, const Real sie,
                                            Real *lambda = nullptr) const {
    const Real temp =
        t_.TemperatureFromDensityInternalEnergy(rho * rho_unit_, sie * sie_unit_, lambda);
    return inv_temp_unit_ * temp;
  }
  PORTABLE_FUNCTION
  Real InternalEnergyFromDensityTemperature(const Real rho, const Real temperature,
                                            Real *lambda = nullptr) const {
    const Real sie = t_.InternalEnergyFromDensityTemperature(
        rho * rho_unit_, temperature * temp_unit_, lambda);
    return inv_sie_unit_ * sie;
  }
  PORTABLE_FUNCTION
  Real PressureFromDensityInternalEnergy(const Real rho, const Real sie,
                                         Real *lambda = nullptr) const {
    const Real P =
        t_.PressureFromDensityInternalEnergy(rho * rho_unit_, sie * sie_unit_, lambda);
    return inv_press_unit_ * P;
  }
  PORTABLE_FUNCTION
  Real SpecificHeatFromDensityInternalEnergy(const Real rho, const Real sie,
                                             Real *lambda = nullptr) const {
    const Real cv = t_.SpecificHeatFromDensityInternalEnergy(rho * rho_unit_,
                                                             sie * sie_unit_, lambda);
    return inv_cv_unit_ * cv;
  }
  PORTABLE_FUNCTION
  Real BulkModulusFromDensityInternalEnergy(const Real rho, const Real sie,
                                            Real *lambda = nullptr) const {
    const Real bmod =
        t_.BulkModulusFromDensityInternalEnergy(rho * rho_unit_, sie * sie_unit_, lambda);
    return inv_bmod_unit_ * bmod;
  }
  PORTABLE_FUNCTION
  Real GruneisenParamFromDensityInternalEnergy(const Real rho, const Real sie,
                                               Real *lambda = nullptr) const {
    const Real gm1 = t_.GruneisenParamFromDensityInternalEnergy(rho * rho_unit_,
                                                                sie * sie_unit_, lambda);
    return gm1;
  }
  PORTABLE_FUNCTION
  Real PressureFromDensityTemperature(const Real rho, const Real temp,
                                      Real *lambda = nullptr) const {
    const Real P =
        t_.PressureFromDensityTemperature(rho * rho_unit_, temp * temp_unit_, lambda);
    return inv_press_unit_ * P;
  }
  PORTABLE_FUNCTION
  Real SpecificHeatFromDensityTemperature(const Real rho, const Real temp,
                                          Real *lambda = nullptr) const {
    const Real cv =
        t_.SpecificHeatFromDensityTemperature(rho * rho_unit_, temp * temp_unit_, lambda);
    return inv_cv_unit_ * cv;
  }
  PORTABLE_FUNCTION
  Real BulkModulusFromDensityTemperature(const Real rho, const Real temp,
                                         Real *lambda = nullptr) const {
    const Real bmod =
        t_.BulkModulusFromDensityTemperature(rho * rho_unit_, temp * temp_unit_, lambda);
    return inv_bmod_unit_ * bmod;
  }
  PORTABLE_FUNCTION
  Real GruneisenParamFromDensityTemperature(const Real rho, const Real temp,
                                            Real *lambda = nullptr) const {
    const Real gm1 = t_.GruneisenParamFromDensityTemperature(rho * rho_unit_,
                                                             temp * temp_unit_, lambda);
    return gm1;
  }

  PORTABLE_FUNCTION
  void DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                            Real *lambda, Real &rho, Real &sie) const {
    t_.DensityEnergyFromPressureTemperature(press * press_unit_, temp * temp_unit_,
                                            lambda, rho, sie);
    rho *= inv_rho_unit_;
    sie *= inv_sie_unit_;
  }

  PORTABLE_FUNCTION
  void PTofRE(const Real rho, const Real sie, Real *lambda, Real &press, Real &temp,
              Real &dpdr, Real &dpde, Real &dtdr, Real &dtde) const {
    t_.PTofRE(rho * rho_unit_, sie * sie_unit_, lambda, press, temp, dpdr, dpde, dtdr,
              dtde);
    press *= inv_press_unit_;
    temp *= inv_temp_unit_;
    dpdr *= inv_dpdr_unit_;
    dpde *= inv_dpde_unit_;
    dtdr *= inv_dtdr_unit_;
    dtde *= inv_dtde_unit_;
  }

  PORTABLE_FUNCTION
  void FillEos(Real &rho, Real &temp, Real &energy, Real &press, Real &cv, Real &bmod,
               const unsigned long output, Real *lambda = nullptr) const {
    // TODO(JMM): Is this general enough? Do I need more switches/scales?
    Real srho = rho_unit_ * rho;
    switch (t_.PreferredInput()) {
    case thermalqs::density | thermalqs::temperature:
      Real sT = temp_unit_ * temp;
      t_.FillEos(srho, sT, energy, press, cv, bmod, output, lambda);
      energy *= inv_sie_unit_;
      break;
    case thermalqs::density | thermalqs::specific_internal_energy:
      Real ssie = sie_unit_ * energy;
      t_.FillEos(srho, temp, ssie, press, cv, bmod, output, lambda);
      break;
    default:
      EOS_ERROR("Didn't find a valid input for ScaledEOS::FillEOS\n");
    }
    press *= inv_press_unit_;
    cv *= inv_cv_unit_;
    bmod *= inv_bmod_unit_;
  }

  PORTABLE_FUNCTION
  void ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                              Real &bmod, Real &dpde, Real &dvdt,
                              Real *lambda = nullptr) const {
    t_.ValuesAtReferenceState(rho, temp, sie, press, cv, bmod, dpde, dvdt, lambda);
    rho *= inv_rho_unit_;
    temp *= inv_temp_unit_;
    sie *= inv_sie_unit_;
    press *= inv_press_unit_;
    cv *= inv_cv_unit_;
    bmod *= inv_bmod_unit_;
    dpde *= inv_dpde_unit_;
    dvdt *= inv_dvdt_unit_;
  }

  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return t_.nlambda(); }

  PORTABLE_FUNCTION
  unsigned long PreferredInput() const { return t_.PreferredInput(); }

  PORTABLE_FUNCTION void PrintParams() const {
    t_.PrintParams();
    printf("Units = %e %e %e %e\n", rho_unit_, sie_unit_, temp_unit_, press_unit_);
  }

 private:
  T t_;

  // JMM: The compiler deletes GetOnDevice if I make these const. So... they're not
  Real rho_unit_, sie_unit_, temp_unit_, press_unit_;
  Real inv_rho_unit_, inv_sie_unit_, inv_temp_unit_, inv_press_unit_;
  Real inv_dpde_unit_, inv_dvdt_unit_, inv_dpdr_unit_, inv_dtdr_unit_, inv_dtde_unit_;
  Real inv_cv_unit_, inv_bmod_unit_;
};

} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_EOS_UNITSYSTEM_HPP_

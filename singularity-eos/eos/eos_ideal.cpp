//------------------------------------------------------------------------------
// © 2021 - 2022. Triad National Security, LLC. All rights reserved.  This
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

#include <singularity-eos/eos/eos.hpp>

#define MYMAX(a, b) a > b ? a : b

namespace singularity {

//---------------------
// Ideal Gas EOS
//---------------------
PORTABLE_FUNCTION Real IdealGas::InternalEnergyFromDensityTemperature(
    const Real rho, const Real temp, Real *lambda) const {
  return MYMAX(0.0, _Cv * temp);
}
PORTABLE_FUNCTION Real IdealGas::TemperatureFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  return MYMAX(0.0, sie / _Cv);
}
PORTABLE_FUNCTION Real IdealGas::PressureFromDensityInternalEnergy(const Real rho,
                                                                   const Real sie,
                                                                   Real *lambda) const {
  return MYMAX(0.0, _gm1 * rho * sie);
}
PORTABLE_FUNCTION Real IdealGas::SpecificHeatFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  return _Cv;
}
PORTABLE_FUNCTION Real IdealGas::BulkModulusFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  return MYMAX(0.0, (_gm1 + 1) * _gm1 * rho * sie);
}

PORTABLE_FUNCTION
Real IdealGas::GruneisenParamFromDensityInternalEnergy(const Real rho, const Real sie,
                                                       Real *lambda) const {
  return _gm1;
}

PORTABLE_FUNCTION void
IdealGas::DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                               Real *lambda, Real &rho, Real &sie) const {
  sie = MYMAX(0.0, _Cv * temp);
  rho = MYMAX(0.0, press / (_gm1 * sie));
}

PORTABLE_FUNCTION
void IdealGas::FillEos(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                       Real &bmod, const unsigned long output, Real *lambda) const {
  if (output & thermalqs::density && output & thermalqs::specific_internal_energy) {
    if (output & thermalqs::pressure || output & thermalqs::temperature) {
      UNDEFINED_ERROR;
    }
    DensityEnergyFromPressureTemperature(press, temp, lambda, rho, sie);
  }
  if (output & thermalqs::pressure && output & thermalqs::specific_internal_energy) {
    if (output & thermalqs::density || output & thermalqs::temperature) {
      UNDEFINED_ERROR;
    }
    sie = InternalEnergyFromDensityTemperature(rho, temp, lambda);
  }
  if (output & thermalqs::temperature && output & thermalqs::specific_internal_energy) {
    sie = press / (_gm1 * rho);
  }
  if (output & thermalqs::pressure) press = PressureFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::temperature)
    temp = TemperatureFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::bulk_modulus)
    bmod = BulkModulusFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::specific_heat)
    cv = SpecificHeatFromDensityInternalEnergy(rho, sie);
}

PORTABLE_FUNCTION
void IdealGas::ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press,
                                      Real &cv, Real &bmod, Real &dpde, Real &dvdt,
                                      Real *lambda) const {
  // use STP: 1 atmosphere, room temperature
  rho = _rho0;
  temp = _T0;
  sie = _sie0;
  press = _P0;
  cv = _Cv;
  bmod = _bmod0;
  dpde = _dpde0;
  dvdt = _dvdt0;
}

PORTABLE_FUNCTION
Real IdealGas::PressureFromDensityTemperature(const Real rho, const Real temp,
                                              Real *lambda) const {
  return MYMAX(0.0, _gm1 * rho * _Cv * temp);
}

PORTABLE_FUNCTION
Real IdealGas::SpecificHeatFromDensityTemperature(const Real rho, const Real temp,
                                                  Real *lambda) const {
  return _Cv;
}
PORTABLE_FUNCTION
Real IdealGas::BulkModulusFromDensityTemperature(const Real rho, const Real temp,
                                                 Real *lambda) const {
  return MYMAX(0.0, (_gm1 + 1) * _gm1 * rho * _Cv * temp);
}

PORTABLE_FUNCTION
Real IdealGas::GruneisenParamFromDensityTemperature(const Real rho, const Real temp,
                                                    Real *lambda) const {
  return _gm1;
}

} // namespace singularity

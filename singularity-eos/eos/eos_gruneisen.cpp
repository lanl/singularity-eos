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

#include <singularity-eos/eos/eos.hpp>

namespace singularity {

PORTABLE_INLINE_FUNCTION Real square(const Real x) { return x * x; }
PORTABLE_INLINE_FUNCTION Real cube(const Real x) { return x * x * x; }

PORTABLE_INLINE_FUNCTION Real Gruneisen::Gamma(const Real rho) const {
  return rho < _rho0 ? _G0 : _G0 * _rho0 / rho + _b * (1 - _rho0 / rho);
}
PORTABLE_INLINE_FUNCTION Real Gruneisen::dPres_drho_e(const Real rho, const Real sie) const {
  if (rho < _rho0) {
    return square(_C0) + Gamma(rho) * sie;
  } else {
    const Real eta = 1 - _rho0 / rho;
    const Real s = _s1 * eta + _s2 * square(eta) + _s3 * cube(eta);
    const Real ds = _s1 + 2 * _s2 * eta + 3 * _s3 * square(eta);
    const Real deta = _rho0 / square(rho);
    const Real dGam = (_b - _G0) * deta;
    const Real P_H = _P0 + square(_C0) * _rho0 * eta / square(1 - s);
    const Real dP_H = square(_C0) * _rho0 / square(1 - s) * deta * (1 + 2 * eta * ds /
      (1 - s));
    const Real E_H = (P_H + _P0) * eta / _rho0 / 2.;
    const Real dE_H = deta * (P_H + _P0) / _rho0 / 2. + eta / _rho0 / 2 * dP_H;
    return dP_H + Gamma(rho) * (sie - E_H) + rho * dGam * (sie - E_H) - 
      rho * Gamma(rho) * dE_H;
  }
}
PORTABLE_FUNCTION Real Gruneisen::InternalEnergyFromDensityTemperature(
    const Real rho, const Real temp, Real *lambda) const {
  return _Cv * (temp - _T0);
}
PORTABLE_FUNCTION Real Gruneisen::TemperatureFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  return _T0 + sie / _Cv;
}
PORTABLE_FUNCTION Real Gruneisen::PressureFromDensityInternalEnergy(
    const Real rho, const Real sie,
    Real *lambda) const { // Will a compiler clean this up for me?
  Real P_H;
  Real E_H;
  if (rho >= _rho0) {
    const Real eta = 1 - _rho0 / rho;
    const Real s = _s1 * eta + _s2 * square(eta) + _s3 * cube(eta);
    P_H = _P0 + square(_C0) * _rho0 * eta / square(1 - s);
    E_H = (P_H + _P0) * eta / _rho0 / 2.;
  } else {
    // This isn't thermodynamically consistent but it's widely used for expansion
    P_H = _P0 + square(_C0) * (rho - _rho0);
    E_H = 0.;
  }
  return P_H + Gamma(rho) * rho * (sie - E_H);
}
PORTABLE_FUNCTION Real Gruneisen::SpecificHeatFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  return _Cv;
}
PORTABLE_FUNCTION Real Gruneisen::BulkModulusFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  // The if statement exists here to avoid the divide by zero
  if (rho < _rho0) {
    return rho * square(_C0) +
           _G0 * (rho * sie + PressureFromDensityInternalEnergy(rho, sie));
  } else {
    const Real dPdr_e = dPres_drho_e(rho, sie);
    const Real dPde_r = rho * Gamma(rho);
    // Thermodynamic identity
    return rho * dPdr_e + PressureFromDensityInternalEnergy(rho, sie) / rho * dPde_r;
  }
}
PORTABLE_FUNCTION
Real Gruneisen::GruneisenParamFromDensityInternalEnergy(const Real rho, const Real sie,
                                                        Real *lambda) const {
  return Gamma(rho);
}
// Below are "unimplemented" routines
PORTABLE_FUNCTION Real Gruneisen::PressureFromDensityTemperature(const Real rho,
                                                                 const Real temp,
                                                                 Real *lambda) const {
  return PressureFromDensityInternalEnergy(
      rho, InternalEnergyFromDensityTemperature(rho, temp));
}
PORTABLE_FUNCTION Real Gruneisen::SpecificHeatFromDensityTemperature(const Real rho,
                                                                     const Real temp,
                                                                     Real *lambda) const {
  return SpecificHeatFromDensityInternalEnergy(
      rho, InternalEnergyFromDensityTemperature(rho, temp));
}
PORTABLE_FUNCTION Real Gruneisen::BulkModulusFromDensityTemperature(const Real rho,
                                                                    const Real temp,
                                                                    Real *lambda) const {
  return BulkModulusFromDensityInternalEnergy(
      rho, InternalEnergyFromDensityTemperature(rho, temp));
}
PORTABLE_FUNCTION
Real Gruneisen::GruneisenParamFromDensityTemperature(const Real rho, const Real temp,
                                                     Real *lambda) const {
  return Gamma(rho);
}
PORTABLE_FUNCTION void Gruneisen::DensityEnergyFromPressureTemperature(
    const Real press, const Real temp, Real *lambda, Real &rho, Real &sie) const {
  sie = _Cv * (temp - _T0);
  // We have a branch at rho0, so we need to decide, based on our pressure, whether we
  // should be above or below rho0
  Real Pref = PressureFromDensityInternalEnergy(_rho0, sie);
  if (press < Pref) {
    rho = (press - _P0 + _C0 * _C0 * _rho0) / (_C0 * _C0 + _G0 * sie);
  } else { // We are in compression; iterate
    auto residual = [&](const Real r) {
      return press - PressureFromDensityInternalEnergy(r, sie);
    };
    Real rho1 = _rho0, res1 = residual(rho1), slope = _G0 * sie + _C0 * _C0, rho2, res2,
         rhom, resm;
    rho2 = (rho > rho1 + 1e-3) ? rho : rho1 + res1 / slope;
    res2 = residual(rho2);
    for (int i = 0; i < 20; ++i) {
      slope = (rho2 - rho1) / (res2 - res1 + 1.0e-10);
      rhom = rho2 - res2 * slope;
      resm = residual(rhom);
      if (resm / press < 1e-8) break;
      rho1 = rho2;
      res1 = res2;
      rho2 = rhom;
      res2 = resm;
    }
    rho = rhom;
  }
}
PORTABLE_FUNCTION void Gruneisen::FillEos(Real &rho, Real &temp, Real &sie, Real &press,
                                          Real &cv, Real &bmod,
                                          const unsigned long output,
                                          Real *lambda) const {
  // The following could be sped up with work!
  if (output & thermalqs::pressure) press = PressureFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::temperature)
    temp = TemperatureFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::bulk_modulus)
    bmod = BulkModulusFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::specific_heat)
    cv = SpecificHeatFromDensityInternalEnergy(rho, sie);
}

// TODO(JMM): pre-cache these rather than recomputing them each time
PORTABLE_FUNCTION
void Gruneisen::ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press,
                                       Real &cv, Real &bmod, Real &dpde, Real &dvdt,
                                       Real *lambda) const {
  rho = _rho0;
  temp = _T0;
  sie = 0;
  press = PressureFromDensityInternalEnergy(rho, sie, lambda);
  cv = SpecificHeatFromDensityInternalEnergy(rho, sie, lambda);
  bmod = BulkModulusFromDensityInternalEnergy(rho, sie, lambda);
  dpde = PressureFromDensityInternalEnergy(rho, sie, lambda) / sie;
  // TODO: chad please fix this one for me. This is wrong.
  Real gm1 = GruneisenParamFromDensityInternalEnergy(rho, sie, lambda) * rho;
  dvdt = gm1 * cv / bmod;
}
} // namespace singularity

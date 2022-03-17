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

#include <singularity-eos/eos/eos.hpp>

namespace singularity {

PORTABLE_INLINE_FUNCTION Real square(const Real x) { return x * x; }
PORTABLE_INLINE_FUNCTION Real cube(const Real x) { return x * x * x; }
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
  const Real top2 = 1.0 - 0.5 * _G0;
  const Real top3 = -0.5 * _b;
  const Real eta = rho / _rho0;
  const Real mu = eta - 1.0;
  Real scale = 1.0;
  const Real G = Gamma(rho);
  if (mu > 0.0) {
    const Real top = 1.0 + mu * (top2 + mu * top3);
    const Real z = mu / eta;
    const Real bot1 = (1.0 - _s1) * mu;
    const Real bot2 = -_s2 * mu * z;
    const Real bot3 = -_s3 * mu * z * z;
    const Real bot = 1.0 + bot1 + bot2 + bot3;
    scale = top / square(bot);
    // if(sie>0.0) G = _G0+_b*mu;
  }
  const Real rho_min = rho < _rho0 ? rho : _rho0;
  return scale * _rho0 * _C0 * _C0 * mu + G * rho_min * sie;
}
PORTABLE_FUNCTION Real Gruneisen::SpecificHeatFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  return _Cv;
}
PORTABLE_FUNCTION Real Gruneisen::BulkModulusFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  const Real mu = rho / _rho0 - 1;
  if (mu <= 0) {
    return rho * square(_C0) +
           _G0 * (rho * sie + PressureFromDensityInternalEnergy(rho, sie));
  } else {
    const Real x = 1 - _rho0 / rho;
    const Real s = _s1 + _s2 * x + _s3 * square(x);
    const Real ds = _s2 + 2 * _s3 * x;
    const Real Pr = _rho0 * square(_C0) * x / square(1 - x * s);
    const Real dPr = -square(_rho0 * _C0) * (1 + x * (s + 2 * x * ds)) / cube(1 - s * x);
    const Real G = _G0 + _b * mu;
    return -dPr / rho * (1 - 0.5 * G * x) - 0.5 * G * Pr +
           G * PressureFromDensityInternalEnergy(rho, sie) / (mu + 1) + rho * _b * sie;
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
    rho = (press + _C0 * _C0 * _rho0) / (_C0 * _C0 + _G0 * sie);
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

#if 0
PORTABLE_FUNCTION void Gruneisen::PTofRE(const Real rho, const Real sie, Real * lambda, Real& press, Real& temp, Real & dpdr, Real & dpde, Real & dtdr, Real & dtde) const
{
    press = PressureFromDensityInternalEnergy(rho,sie);
    temp = TemperatureFromDensityInternalEnergy(rho,sie);
    const Real u = rho/_rho0 - 1.0;
    if(u>0) {
        const Real b0 = 1.0-_s1;
	const Real b3 = _s2+_s3;
  const Real d = 1.0+(b0+1)
    dpdr = (fu/square(d)-2.0*aa*du/d)/_rho0+_b*sie;
	dpde = (_G0+_b*u)*_rho0;
    } else {
      dpdr = _C0*_C0;
      dpde = _G0*_rho0;
    }
}
#endif

} // namespace singularity

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

#include <cmath>
#include <singularity-eos/eos/eos.hpp>

namespace singularity {

PORTABLE_INLINE_FUNCTION Real square(const Real x) { return x * x; }
PORTABLE_INLINE_FUNCTION Real cube(const Real x) { return x * x * x; }
PORTABLE_INLINE_FUNCTION Real fourth(const Real x) { return x * x * x * x; }
PORTABLE_INLINE_FUNCTION Real fifth(const Real x) { return x * x * x * x * x; }

constexpr Real onethird = 1.0 / 3.0;

PORTABLE_INLINE_FUNCTION Real DavisReactants::Ps(const Real rho) const {
  const Real y = 1.0 - _rho0 / rho;
  const Real phat = 0.25 * _A * _A / _B * _rho0;
  const Real b4y = 4.0 * _B * y;

  if (rho >= _rho0) {
    return phat *
           (b4y + 0.5 * (square(b4y) + onethird * (cube(b4y) + _C * fourth(b4y) * 0.25)) +
            square(y) / fourth(1 - y));
  } else {
    return phat * (std::exp(b4y) - 1.0);
  }
}
PORTABLE_INLINE_FUNCTION Real DavisReactants::Es(const Real rho) const {
  const Real y = 1 - _rho0 / rho;
  const Real phat = 0.25 * _A * _A / _B * _rho0;
  const Real b4y = 4 * _B * y;
  Real e_s;
  if (y > 0.0) {
    const Real z = rho / _rho0 - 1;
    e_s = 0.5 * y * b4y *
              (1.0 + onethird * b4y * (1.0 + 0.25 * b4y * (1.0 + _C * 0.2 * b4y))) +
          onethird * cube(z);
  } else {
    e_s = -y - (1.0 - std::exp(b4y)) / (4.0 * _B);
  }
  return _e0 + _P0 * (1.0 / _rho0 - 1.0 / rho) + phat / _rho0 * e_s;
}
PORTABLE_INLINE_FUNCTION Real DavisReactants::Ts(const Real rho) const {
  if (rho >= _rho0) {
    const Real y = 1 - _rho0 / rho;
    return _T0 * std::exp(-_Z * y) * std::pow(_rho0 / rho, -_G0 - _Z);
  } else {
    return _T0 * std::pow(_rho0 / rho, -_G0);
  }
}
PORTABLE_INLINE_FUNCTION Real DavisReactants::Gamma(const Real rho) const {
  if (rho >= _rho0) {
    const Real y = 1 - _rho0 / rho;
    return _G0 + _Z * y;
  } else {
    return _G0;
  }
}
PORTABLE_FUNCTION Real DavisReactants::InternalEnergyFromDensityTemperature(
    const Real rho, const Real temp, Real *lambda) const {
  const Real t_s = Ts(rho);
  return Es(rho) +
         _Cv0 * t_s / (1.0 + _alpha) * (std::pow(temp / t_s, 1.0 + _alpha) - 1.0);
}
PORTABLE_FUNCTION Real DavisReactants::PressureFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  return Ps(rho) + Gamma(rho) * rho * (sie - Es(rho));
}
PORTABLE_FUNCTION Real DavisReactants::TemperatureFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  const Real es = Es(rho);
  const Real tmp = std::pow((1.0 + _alpha) / (Ts(rho) * _Cv0) * (sie - es) + 1.0,
                            1.0 / (1.0 + _alpha));
  if (tmp > 0) return Ts(rho) * tmp;
  return Ts(rho) + (sie - es) / _Cv0; // This branch is a negative temperature
}
PORTABLE_FUNCTION Real DavisReactants::SpecificHeatFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  return _Cv0 / std::pow((1 + _alpha) / (Ts(rho) * _Cv0) * (sie - Es(rho)) + 1,
                         -_alpha / (1 + _alpha));
}
PORTABLE_FUNCTION Real DavisReactants::BulkModulusFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  const Real y = 1 - _rho0 / rho;
  const Real phat = 0.25 * _A * _A / _B * _rho0;
  const Real b4y = 4 * _B * y;
  const Real gamma = Gamma(rho);
  const Real esv = -Ps(rho);
  const Real psv =
      (rho >= _rho0)
          ? -phat * _rho0 *
                (4 * _B * (1 + b4y + 0.5 * (square(b4y) + _C / 3 * cube(b4y))) +
                 3 * y / fourth(1 - y) + 4 * square(y) / fifth(1 - y))
          : -phat * 4 * _B * _rho0 * std::exp(b4y);
  const Real gammav = (rho >= _rho0) ? _Z * _rho0 : 0.0;
  return -(psv + (sie - Es(rho)) * rho * (gammav - gamma * rho) - gamma * rho * esv) /
         rho;
}
PORTABLE_FUNCTION
Real DavisReactants::GruneisenParamFromDensityInternalEnergy(const Real rho,
                                                             const Real sie,
                                                             Real *lambda) const {
  return Gamma(rho);
}
// Below are "unimplemented" routines
PORTABLE_FUNCTION Real DavisReactants::PressureFromDensityTemperature(
    const Real rho, const Real temp, Real *lambda) const {
  return PressureFromDensityInternalEnergy(
      rho, InternalEnergyFromDensityTemperature(rho, temp));
}
PORTABLE_FUNCTION Real DavisReactants::SpecificHeatFromDensityTemperature(
    const Real rho, const Real temp, Real *lambda) const {
  return SpecificHeatFromDensityInternalEnergy(
      rho, InternalEnergyFromDensityTemperature(rho, temp));
}
PORTABLE_FUNCTION Real DavisReactants::BulkModulusFromDensityTemperature(
    const Real rho, const Real temp, Real *lambda) const {
  return BulkModulusFromDensityInternalEnergy(
      rho, InternalEnergyFromDensityTemperature(rho, temp));
}
PORTABLE_FUNCTION
Real DavisReactants::GruneisenParamFromDensityTemperature(const Real rho, const Real temp,
                                                          Real *lambda) const {
  return Gamma(rho);
}
PORTABLE_FUNCTION void DavisReactants::DensityEnergyFromPressureTemperature(
    const Real press, const Real temp, Real *lambda, Real &rho, Real &sie) const {
  // First, solve P=P(rho,T) for rho.  Note P(rho,e) has an sie-es term, which is only a
  // function of T
  auto residual = [&](const Real r) {
    return press - (Ps(r) + Gamma(r) * r * _Cv0 * Ts(r) / (1 + _alpha) *
                                (std::pow(temp / Ts(r), 1 + _alpha) - 1.0));
  };
  Real rho1 = _rho0, res1 = residual(rho1);
  Real slope =
      _Cv0 * _G0 *
      ((_alpha * _G0 - 1.0) * temp * std::pow(temp / _T0, _alpha) + _T0 + _G0 * _T0) /
      (1 + _alpha);
  Real rho2 = rho1 + res1 / slope, res2 = residual(rho2);
  Real rhom, resm;
  for (int i = 1; i < 20; ++i) {
    slope = (rho2 - rho1) / (res2 - res1);
    rhom = rho2 - res2 * slope;
    resm = residual(rhom);
    if (resm / press < 1e-8) break;
    rho1 = rho2;
    res1 = res2;
    rho2 = rhom;
    res2 = resm;
  }
  rho = rhom;
  sie = InternalEnergyFromDensityTemperature(rho, temp);
}
PORTABLE_FUNCTION
void DavisReactants::FillEos(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                             Real &bmod, const unsigned long output, Real *lambda) const {
  if (output & thermalqs::pressure) press = PressureFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::temperature)
    temp = TemperatureFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::bulk_modulus)
    bmod = BulkModulusFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::specific_heat)
    cv = SpecificHeatFromDensityInternalEnergy(rho, sie);
}
// TODO: Chad please decide if this is sane
// TODO(JMM): Pre-cache values instead of computing inline
PORTABLE_FUNCTION
void DavisReactants::ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press,
                                            Real &cv, Real &bmod, Real &dpde, Real &dvdt,
                                            Real *lambda) const {
  rho = _rho0;
  temp = _T0;
  sie = InternalEnergyFromDensityTemperature(_rho0, _T0);
  press = PressureFromDensityInternalEnergy(_rho0, sie);
  cv = SpecificHeatFromDensityInternalEnergy(_rho0, sie);
  bmod = BulkModulusFromDensityInternalEnergy(_rho0, sie);
  dpde = Gamma(_rho0) * _rho0;
  // chad please fix this. it's wrong
  Real gm1 = GruneisenParamFromDensityTemperature(_rho0, _T0) * _rho0;
  dvdt = gm1 * cv / bmod;
}

#if 0
PORTABLE_FUNCTION void DavisReactants::PTofRE(const Real rho, const Real sie, Real * lambda, Real& press, Real& temp, Real & dpdr, Real & dpde, Real & dtdr, Real & dtde) const
{
    press = PressureFromDensityInternalEnergy(rho,sie,lambda);
    temp = TemperatureFromDensityInternalEnergy(rho,sie,lambda);
    const Real y = 1-_rho0/rho;
    const Real phat = 0.25*_A*_A/_B*_rho0;
    const Real b4y = 4*_B*y;
    const Real gamma = Gamma(rho);
    const Real esv = -Ps(rho);
    const Real psv = (rho>=_rho0)? -phat*_rho0*(4*_B*(1+b4y+0.5*(square(b4y)+_C/3*cube(b4y)))+3*y/fourth(1-y)+4*square(y)/fifth(1-y)) : -phat*4*_B*_rho0*std::exp(b4y);
    const Real gammav = (rho>=_rho0)? _Z*_rho0 : 0.0;
    dpdr = -(psv+(sie-Es(rho))*rho*(gammav-gamma*rho)-gamma*rho*esv)/square(rho);
    dpde = gamma*rho;
    const Real tmp = (1.0+_alpha)/(_Cv0*Ts(rho))*(sie-Es(rho)) + 1.0;
    const Real Tsv = -Ts(rho)*rho*(_G0+(rho>=_rho0)?_Z*y:0.0);
    const Real Tv = std::pow(tmp,1.0/(1.0+_alpha))*(Tsv+Tsv/Ts(rho)*(Es(rho)-sie)-esv)/tmp/_Cv0;
    
    dtdr = -Tv/square(rho);
    dtde = 1.0/SpecificHeatFromDensityInternalEnergy(rho,sie);
}
#endif

PORTABLE_INLINE_FUNCTION Real DavisProducts::F(const Real rho) const {
  const Real vvc = 1.0 / (rho * _vc);
  return 2.0 * _a / (std::pow(vvc, 2 * _n) + 1.0);
}
PORTABLE_INLINE_FUNCTION Real DavisProducts::Es(const Real rho) const {
  const Real vvc = 1 / (rho * _vc);
  const Real ec = _pc * _vc / (_k - 1.0 + _a);
  // const Real de = ecj-(Es(rho0)-_E0);
  return ec * std::pow(0.5 * (std::pow(vvc, _n) + std::pow(vvc, -_n)), _a / _n) /
         std::pow(vvc, _k - 1.0 + _a);
}
PORTABLE_INLINE_FUNCTION Real DavisProducts::Ps(const Real rho) const {
  const Real vvc = 1 / (rho * _vc);
  return _pc * std::pow(0.5 * (std::pow(vvc, _n) + std::pow(vvc, -_n)), _a / _n) /
         std::pow(vvc, _k + _a) * (_k - 1.0 + F(rho)) / (_k - 1.0 + _a);
}
PORTABLE_INLINE_FUNCTION Real DavisProducts::Ts(const Real rho) const {
  const Real vvc = 1 / (rho * _vc);
  return std::pow(2.0, -_a * _b / _n) * _pc * _vc / (_Cv * (_k - 1 + _a)) *
         std::pow(0.5 * (std::pow(vvc, _n) + std::pow(vvc, -_n)), _a / _n * (1 - _b)) /
         std::pow(vvc, _k - 1.0 + _a * (1 - _b));
}
PORTABLE_INLINE_FUNCTION Real DavisProducts::Gamma(const Real rho) const {
  return _k - 1.0 + (1.0 - _b) * F(rho);
}
PORTABLE_FUNCTION Real DavisProducts::InternalEnergyFromDensityTemperature(
    const Real rho, const Real temp, Real *lambda) const {
  return _Cv * (temp - Ts(rho)) + Es(rho);
}
PORTABLE_FUNCTION Real DavisProducts::PressureFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  return Ps(rho) + rho * Gamma(rho) * (sie - Es(rho));
}
PORTABLE_FUNCTION Real DavisProducts::TemperatureFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  return Ts(rho) + (sie - Es(rho)) / _Cv;
}
PORTABLE_FUNCTION Real DavisProducts::SpecificHeatFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  return _Cv;
}
PORTABLE_FUNCTION Real DavisProducts::BulkModulusFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  const Real vvc = 1 / (rho * _vc);
  const Real Fx = -4 * _a * std::pow(vvc, 2 * _n - 1) / square(1 + std::pow(vvc, 2 * _n));
  const Real tmp = std::pow(0.5 * (std::pow(vvc, _n) + std::pow(vvc, -_n)), _a / _n) /
                   std::pow(vvc, _k + _a);
  const Real tmp_x =
      0.5 * _a * (std::pow(vvc, _n - 1) - std::pow(vvc, -_n - 1)) *
          std::pow(0.5 * (std::pow(vvc, _n) + std::pow(vvc, -_n)), _a / _n - 1) /
          std::pow(vvc, _k + _a) -
      (_k + _a) * tmp / vvc;
  const Real psv = _pc / (_k - 1 + _a) * (tmp * Fx + (_k - 1 + F(rho)) * tmp_x) / _vc;
  // const Real esv = _pc*_vc/(_k-1+_a)*(tmp+vvc*tmp_x)/_vc;
  const Real esv = _pc / (_k - 1 + _a) * (tmp + vvc * tmp_x);
  const Real gamma = Gamma(rho);
  const Real gammav = (1 - _b) * Fx * _vc;
  return -(psv + (sie - Es(rho)) * rho * (gammav - gamma * rho) - gamma * rho * esv) /
         rho;
}
PORTABLE_FUNCTION
Real DavisProducts::GruneisenParamFromDensityInternalEnergy(const Real rho,
                                                            const Real sie,
                                                            Real *lambda) const {
  return Gamma(rho);
}
// Below are "unimplemented" routines
PORTABLE_FUNCTION Real DavisProducts::PressureFromDensityTemperature(const Real rho,
                                                                     const Real temp,
                                                                     Real *lambda) const {
  return PressureFromDensityInternalEnergy(
      rho, InternalEnergyFromDensityTemperature(rho, temp));
}
PORTABLE_FUNCTION Real DavisProducts::SpecificHeatFromDensityTemperature(
    const Real rho, const Real temp, Real *lambda) const {
  return SpecificHeatFromDensityInternalEnergy(
      rho, InternalEnergyFromDensityTemperature(rho, temp));
}
PORTABLE_FUNCTION Real DavisProducts::BulkModulusFromDensityTemperature(
    const Real rho, const Real temp, Real *lambda) const {
  return BulkModulusFromDensityInternalEnergy(
      rho, InternalEnergyFromDensityTemperature(rho, temp));
}
PORTABLE_FUNCTION
Real DavisProducts::GruneisenParamFromDensityTemperature(const Real rho, const Real temp,
                                                         Real *lambda) const {
  return Gamma(rho);
}
PORTABLE_FUNCTION void DavisProducts::DensityEnergyFromPressureTemperature(
    const Real press, const Real temp, Real *lambda, Real &rho, Real &sie) const {
  auto residual = [&](const Real r) {
    return press - (Ps(r) + Gamma(r) * r * _Cv * (temp - Ts(r)));
  };
  Real rho1 = 1.0 / _vc, res1 = residual(rho1);
  Real slope;
  Real rho2 = rho1 * 2.0, res2 = residual(rho2);
  Real rhom, resm;
  for (int i = 1; i < 20; ++i) {
    slope = (rho2 - rho1) / (res2 - res1);
    rhom = rho2 - res2 * slope;
    resm = residual(rhom);
    if (resm / press < 1e-8) break;
    rho1 = rho2;
    res1 = res2;
    rho2 = rhom;
    res2 = resm;
  }
  rho = rhom;
  sie = InternalEnergyFromDensityTemperature(rho, temp);
}
PORTABLE_FUNCTION
void DavisProducts::FillEos(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                            Real &bmod, const unsigned long output, Real *lambda) const {
  if (output & thermalqs::pressure) press = PressureFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::temperature)
    temp = TemperatureFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::bulk_modulus)
    bmod = BulkModulusFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::specific_heat)
    cv = SpecificHeatFromDensityInternalEnergy(rho, sie);
}
// TODO: pre-cache values instead of computing them
// TODO: chad please decide if these choices are sane
PORTABLE_FUNCTION
void DavisProducts::ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press,
                                           Real &cv, Real &bmod, Real &dpde, Real &dvdt,
                                           Real *lambda) const {
  rho = 1.0 / _vc;
  sie = 2. * Es(rho);
  temp = TemperatureFromDensityInternalEnergy(rho, sie);
  press = PressureFromDensityInternalEnergy(rho, sie);
  cv = SpecificHeatFromDensityInternalEnergy(rho, sie);
  bmod = BulkModulusFromDensityInternalEnergy(rho, sie);
  dpde = Gamma(rho) * rho;
  // chad please fix this. it's wrong
  Real gm1 = GruneisenParamFromDensityTemperature(rho, temp) * rho;
  dvdt = gm1 * cv / bmod;
}
#if 0
PORTABLE_FUNCTION void DavisProducts::PTofRE(const Real rho, const Real sie, Real * lambda, Real& press, Real& temp, Real & dpdr, Real & dpde, Real & dtdr, Real & dtde) const
{
    press = PressureFromDensityInternalEnergy(rho,sie,lambda);
    temp = TemperatureFromDensityInternalEnergy(rho,sie,lambda);

    const Real vvc = 1/(rho*_vc);
    const Real Fx = -4*_a*std::pow(vvc,2*_n-1)/square(1+std::pow(vvc,2*_n));
    const Real tmp = std::pow(0.5*(std::pow(vvc,_n)+std::pow(vvc,-_n)),_a/_n)/std::pow(vvc,_k+_a);
    const Real tmp_x = 0.5*_a*(std::pow(vvc,_n-1)-std::pow(vvc,-_n-1))*std::pow(0.5*(std::pow(vvc,_n)+std::pow(vvc,-_n)),_a/_n-1)/std::pow(vvc,_k+_a)-(_k+_a)*tmp/vvc;
    const Real psv = _pc/(_k-1+_a)*(tmp*Fx+(_k-1+F(rho))*tmp_x)/_vc;
    //const Real esv = _pc*_vc/(_k-1+_a)*(tmp+vvc*tmp_x)/_vc;
    const Real esv = _pc/(_k-1+_a)*(tmp+vvc*tmp_x);
    const Real gamma = Gamma(rho);
    const Real gammav = (1-_b)*Fx*_vc;
    dpde = gamma*rho;
    dpdr = (psv+(sie-Es(rho))*rho*(gammav-gamma*rho)-gamma*rho*esv)/square(rho);
    dtde = 1.0/_Cv;
    const Real tsv = Ts(rho)*rho*(_a*(1.0-_b)*(std::pow(vvc,2.0*_n)-1.0)/(std::pow(vvc,2.0*_n)+1.0)-(_k-1.0+_a*(1.0-_b)));
    dtdr = -(tsv-esv/_Cv)/square(rho);
}
#endif

} // namespace singularity

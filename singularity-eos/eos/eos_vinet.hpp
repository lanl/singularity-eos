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

#ifndef _SINGULARITY_EOS_EOS_VINET_HPP_
#define _SINGULARITY_EOS_EOS_VINET_HPP_

#include <cmath>
#include <cstdio>

#include <limits>

#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/math_utils.hpp>
#include <singularity-eos/base/root-finding-1d/root_finding.hpp>
#include <singularity-eos/eos/eos_base.hpp>

namespace singularity {

using namespace eos_base;

// COMMENT: This is a complete, thermodynamically consistent Vinet
// EOS created at Sandia National Laboratories by AEM.
// J. Appl. Phys. 119, 015904 (2016)
// It is particularly suited for high-pressure applications.
class Vinet : public EosBase<Vinet> {
 public:
  Vinet() = default;
  // Constructor 
  PORTABLE_INLINE_FUNCTION
  Vinet(const Real rho0, const Real T0, const Real B0, const Real BP0, const Real A0, 
	    const Real Cv0, const Real E0, const Real S0, const Real d2, const Real d3)
      : _rho0(rho0), _T0(T0), _B0(B0), _BP0(BP0), _A0(A0), _Cv0(Cv0), _E0(E0), _S0(S0), 
	_d2(d2), _d3(d3) {
      CheckVinet();
  }

  Vinet GetOnDevice() { return *this; }
  PORTABLE_INLINE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temp, Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityTemperature(
      const Real rho, const Real temp, Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityTemperature(
      const Real rho, const Real temp, Real *lambda = nullptr) const {
    return _Cv0;
  }
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const {
    return _Cv0;
  }
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityTemperature(
      const Real rho, const Real temp, Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temp, Real *lambda = nullptr) const {
    return _A0 * _B0 / _Cv0 / rho;
  }
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const {
    return _A0 * _B0 / _Cv0 / rho ;
  }
  PORTABLE_INLINE_FUNCTION void FillEos(Real &rho, Real &temp, Real &energy, Real &press,
                                        Real &cv, Real &bmod, const unsigned long output,
                                        Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION
  void ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                              Real &bmod, Real &dpde, Real &dvdt,
                              Real *lambda = nullptr) const;
  // Generic functions provided by the base class. These contain e.g. the vector
  // overloads that use the scalar versions declared here
  SG_ADD_BASE_CLASS_USINGS(Vinet)
  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return 0; }
  static constexpr unsigned long PreferredInput() { return _preferred_input; }
  PORTABLE_INLINE_FUNCTION void PrintParams() const {
    static constexpr char st[]{"Vinet Params: "};
    printf("%s rho0:%e T0:%e B0:%e BP0:%e\n  A0:%e Cv0:%e E0:%e S0:%e\n  d2:%e d3:%e \n",
           st, _rho0, _T0, _B0, _BP0, _A0, _Cv0, _E0, _S0, _d2, _d3);
  }
  PORTABLE_INLINE_FUNCTION void
  DensityEnergyFromPressureTemperature(const Real press, const Real temp, Real *lambda,
                                       Real &rho, Real &sie) const;
  inline void Finalize() {}
  static std::string EosType() { return std::string("Vinet"); }
  static std::string EosPyType() { return EosType(); }

 private:
  Real _rho0, _T0, _B0, _BP0, _A0, _Cv0, _E0, _S0, _d2, _d3;
  // static constexpr const char _eos_type[] = {"Vinet"};
  Real *_VIP;
  PORTABLE_INLINE_FUNCTION
  Real dPres_drho_e(const Real rho, const Real sie) const;
  static constexpr const unsigned long _preferred_input =
      thermalqs::density | thermalqs::specific_internal_energy;
  PORTABLE_INLINE_FUNCTION void CheckVinet() const;
  PORTABLE_INLINE_FUNCTION void Vinet_F_DT_func(const Real rho, const Real T, 
		  				Real *output) const;
};

PORTABLE_INLINE_FUNCTION Real Vinet::dPres_drho_e(const Real rho_in,
                                                  const Real sie) const {
  using namespace math_utils;
  const Real rho = std::min(rho_in, _rho_max);
  if (rho < _rho0) {
    return pow<2>(_C0) + Gamma(rho) * sie;
  } else {
    const Real eta = 1 - _rho0 / rho;
    const Real s = _s1 * eta + _s2 * pow<2>(eta) + _s3 * pow<3>(eta);
    const Real ds = _s1 + 2 * _s2 * eta + 3 * _s3 * pow<2>(eta);
    const Real deta = _rho0 / pow<2>(rho);
    const Real dGam = (_b - _G0) * deta;
    const Real P_H = _P0 + pow<2>(_C0) * _rho0 * eta / pow<2>(1 - s);
    const Real dP_H = math_utils::pow<2>(_C0) * _rho0 / pow<2>(1 - s) * deta *
                      (1 + 2 * eta * ds / (1 - s));
    const Real E_H = (P_H + _P0) * eta / _rho0 / 2.;
    const Real dE_H = deta * (P_H + _P0) / _rho0 / 2. + eta / _rho0 / 2 * dP_H;
    return dP_H + Gamma(rho) * (sie - E_H) + rho * dGam * (sie - E_H) -
           rho * Gamma(rho) * dE_H;
  }
}

PORTABLE_INLINE_FUNCTION void Vinet::CheckVinet() const {

  Real *VIP;
  int ind;

  if (VP[VR0] < (EOS_REAL) 0.0){
    printf("testing VP[VR0]\n");
    EOS_ERROR("Required Vinet model parameter rho0 not set");
    return;
  }
  if (VP[VT0] < (EOS_REAL) 0.0){
    printf("testing VP[VT0]\n");
    EOS_ERROR("Required Vinet model parameter T0 not set");
    return;
  }
  if (VP[VBT0] < (EOS_REAL) 0.0){
    printf("testing VP[VBT0]\n");
    EOS_ERROR("Required Vinet model parameter B0 not set");
    return;
  }
  if (VP[VBTP0] < 1.0){
    printf("testing VP[VBTP0]\n");
    EOS_ERROR("Required Vinet model parameter BP0 not set or too small (<1)");
    return;
  }
  if (VP[VA0] < 0.0){
    printf("testing VP[VA0]\n");
    EOS_ERROR("Required Vinet model parameter A0 not set");
    return;
  }
  if (VP[VCV0] < 0.0){
    printf("testing VP[VCV0]\n");
    EOS_ERROR("Required Vinet model parameter Cv0 not set");
    return;
  }

  VIP[0] = _B0/_rho0 + pow(_A0*_B0/_rho0,2.0)*_T0/_Cv0;      // sound speed squared
  VIP[1] = 3.0/2.0*(_BP0-1.0); // exponent
  VIP[2] = 1.0; // prefactor f0
  VIP[3] = 0.0; // prefactor f1

  for (ind = 4; ind < 43; ind++) VIP[ind] = VP[VD2+ind-4];
  ind = 42;
  while (VIP[ind] == 0.0) ind--;
  while (ind >= 4) {
    ind--;
    VIP[ind] = VIP[ind] - (ind -1)/VIP[1]*VIP[ind+1];
  }
}


PORTABLE_INLINE_FUNCTION void Vinet::Vinet_F_DT_func(const Real rho, const Real T, 
							Real *output) const {
  int pref0vp = VD2-2, pref0vip = 2, ind;
  Real x, x2inv, onemx, etatimes1mx, expetatimes1mx, temp;
  Real sumP=0.0, sumB=0.0, sumE=0.0;
  Real entropy, energy, pressure, dpdrho, dpdt, dedt, dedrho, soundspeed;

  ind = 40;
  while (VIP[pref0vip+ind] == 0.0) ind--;

  x = pow(_rho0/rho,1.0/3.0);  /*rho-dependent*/
  x2inv = 1.0/pow(x,2.0);
  onemx = 1.0-x;
  etatimes1mx = VIP[1]*onemx;
  expetatimes1mx = exp(etatimes1mx);

  while (ind >= 2){
    sumP = VP[pref0vp+ind]+onemx*sumP;
    sumB = VP[pref0vp+ind]*ind+onemx*sumB;
    sumE = VIP[pref0vip+ind]+onemx*sumE;
    ind--;
  }
  while (ind >= 0){
    sumE = VIP[pref0vip+ind]+onemx*sumE;
    ind--;
  }
  sumP = 1.0+sumP*pow(onemx,2.0);
  sumB = sumB*onemx;
  sumE = VIP[1]*onemx*sumE;

  /* T0 isotherm */
  energy = 9.0*_B0/pow(VIP[1],2.0)/_rho0(VIP[pref0vip]-expetatimes1mx*(VIP[pref0vip]-sumE))-_A0*_B0*_T0*(1.0-_rho0/rho);
  pressure = 3.0*VP[VBT0]*x2inv*onemx*expetatimes1mx*sumP;
  temp = (1.0+onemx*(VIP[1]*x+1.0))*sumP+x*onemx*sumB;
  temp = VP[VBT0]/VP[VR0]*x*expetatimes1mx*temp;

  /* Go to required temperature */
  energy = energy + VP[VCV0]*(T-VP[VT0])+VP[VESFT];
  pressure = pressure + VP[VA0]*VP[VBT0]*(T-VP[VT0]);
  dpdrho = temp;
  dpdt = VP[VA0]*VP[VBT0];
  dedt = VP[VCV0];
  dedrho = (pressure -T*VP[VA0]*VP[VBT0])/pow(rho,2.0);
  if (T < 0.0){
    printf("warning!! negative temperature");
    entropy = 0.0;
  } else
    entropy = VP[VA0]*VP[VBT0]*(1.0/rho-1/VP[VR0])+VP[VCV0]*log(T/VP[VT0])+VP[VSSFT];
  soundspeed = MAX(0.0,dpdrho+T*pow(dpdt/rho,2.0)/VP[VCV0]);
  soundspeed = sqrt(soundspeed);

  output[0] = energy;
  output[1] = pressure;
  output[2] = dpdrho;
  output[3] = dpdt;
  output[4] = dedt;
  output[5] = dedrho;
  output[6] = entropy;
  output[7] = soundspeed;

  return;
}


PORTABLE_INLINE_FUNCTION Real Vinet::PressureFromDensityInternalEnergy(
    const Real rho_in, const Real sie, Real *lambda) const {
  using namespace math_utils;
  const Real rho = std::min(rho_in, _rho_max);
  Real P_H;
  Real E_H;
  if (rho >= _rho0) {
    const Real eta = 1 - _rho0 / rho;
    const Real s = _s1 * eta + _s2 * pow<2>(eta) + _s3 * pow<3>(eta);
    P_H = _P0 + pow<2>(_C0) * _rho0 * eta / pow<2>(1 - s);
    E_H = (P_H + _P0) * eta / _rho0 / 2.;
  } else {
    // This isn't thermodynamically consistent but it's widely used for expansion
    P_H = _P0 + pow<2>(_C0) * (rho - _rho0);
    E_H = 0.;
  }
  return P_H + Gamma(rho) * rho * (sie - E_H);
}
PORTABLE_INLINE_FUNCTION Real Vinet::BulkModulusFromDensityInternalEnergy(
    const Real rho_in, const Real sie, Real *lambda) const {
  using namespace gruneisen_utils;
  const Real rho = std::min(rho_in, _rho_max);
  // The if statement exists here to avoid the divide by zero
  if (rho < _rho0) {
    return rho * math_utils::pow<2>(_C0) +
           _G0 * (rho * sie + PressureFromDensityInternalEnergy(rho, sie));
  } else {
    const Real dPdr_e = dPres_drho_e(rho, sie);
    const Real dPde_r = rho * Gamma(rho);
    // Thermodynamic identity
    return rho * dPdr_e + PressureFromDensityInternalEnergy(rho, sie) / rho * dPde_r;
  }
}
// Below are "unimplemented" routines
PORTABLE_INLINE_FUNCTION Real Vinet::PressureFromDensityTemperature(
    const Real rho_in, const Real temp, Real *lambda) const {
  const Real rho = std::min(rho_in, _rho_max);
  return PressureFromDensityInternalEnergy(
      rho, InternalEnergyFromDensityTemperature(rho, temp));
}
PORTABLE_INLINE_FUNCTION Real Vinet::BulkModulusFromDensityTemperature(
    const Real rho_in, const Real temp, Real *lambda) const {
  const Real rho = std::min(rho_in, _rho_max);
  return BulkModulusFromDensityInternalEnergy(
      rho, InternalEnergyFromDensityTemperature(rho, temp));
}
PORTABLE_INLINE_FUNCTION void Vinet::DensityEnergyFromPressureTemperature(
    const Real press, const Real temp, Real *lambda, Real &rho, Real &sie) const {
  sie = _Cv * (temp - _T0);
  // We have a branch at rho0, so we need to decide, based on our pressure, whether we
  // should be above or below rho0
  Real Pref = PressureFromDensityInternalEnergy(_rho0, sie);
  if (press < Pref) {
    rho = (press - _P0 + _C0 * _C0 * _rho0) / (_C0 * _C0 + _G0 * sie);
  } else { // We are in compression; iterate
    auto PofRatE = [&](const Real r) {
      return PressureFromDensityInternalEnergy(r, sie);
    };
    using RootFinding1D::regula_falsi;
    using RootFinding1D::Status;
    RootFinding1D::RootCounts counts;
    auto status =
        regula_falsi(PofRatE, press, _rho0, 1.0e-5, 1.0e3, 1.0e-8, 1.0e-8, rho, counts);
    if (status != Status::SUCCESS) {
      // Root finder failed even though the solution was bracketed... this is an error
      EOS_ERROR("Vinet::DensityEnergyFromPressureTemperature: "
                "Root find failed to find a solution given P, T\n");
    }
  }
}
PORTABLE_INLINE_FUNCTION void Vinet::FillEos(Real &rho_in, Real &temp, Real &sie,
                                                 Real &press, Real &cv, Real &bmod,
                                                 const unsigned long output,
                                                 Real *lambda) const {
  // The following could be sped up with work!
  const unsigned long input = ~output;
  if (thermalqs::temperature & input && thermalqs::pressure & input) {
    DensityEnergyFromPressureTemperature(press, temp, lambda, rho_in, sie);
  } else if (thermalqs::density & output ||
             thermalqs::specific_internal_energy & output) {
    // Error out on density or energy output because they're currently required as inputs
    EOS_ERROR("Vinet FillEos: Density and energy are currently required inputs "
              "except when pressure and temperature are inputs.\n");
  }
  const Real rho = std::min(rho_in, _rho_max);
  if (output & thermalqs::pressure) press = PressureFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::temperature)
    temp = TemperatureFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::bulk_modulus)
    bmod = BulkModulusFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::specific_heat)
    cv = SpecificHeatFromDensityInternalEnergy(rho, sie);
}

// TODO(JMM): pre-cache these rather than recomputing them each time
PORTABLE_INLINE_FUNCTION
void Vinet::ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press,
                                       Real &cv, Real &bmod, Real &dpde, Real &dvdt,
                                       Real *lambda) const {
  rho = _rho0;
  temp = _T0;
  sie = _E0;
  press = PressureFromDensityInternalEnergy(rho, sie, lambda);
  cv = SpecificHeatFromDensityInternalEnergy(rho, sie, lambda);
  bmod = BulkModulusFromDensityInternalEnergy(rho, sie, lambda);
  dpde = PressureFromDensityInternalEnergy(rho, sie, lambda) / sie;
  // TODO: chad please fix this one for me. This is wrong.
  Real gm1 = GruneisenParamFromDensityInternalEnergy(rho, sie, lambda) * rho;
  dvdt = gm1 * cv / bmod;
}
} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_VINET_HPP_

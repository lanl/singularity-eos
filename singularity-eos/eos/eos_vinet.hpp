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
	    const Real Cv0, const Real E0, const Real S0, const Real *expconsts)
      : _rho0(rho0), _T0(T0), _B0(B0), _BP0(BP0), _A0(A0), _Cv0(Cv0), _E0(E0), _S0(S0), 
	_expconsts(expconsts) {
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
  PORTABLE_INLINE_FUNCTION Real EntropyFromDensityTemperature(
      const Real rho, const Real temp, Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION Real EntropyFromDensityInternalEnergy(
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
    printf("%s rho0:%e T0:%e B0:%e BP0:%e\n  A0:%e Cv0:%e E0:%e S0:%e\n\n",
           st, _rho0, _T0, _B0, _BP0, _A0, _Cv0, _E0, _S0);
    for( int i=0; i<=38; i++) {
	    if(_d2to40[i] > 0.0) printf("d%i:%e\n",i+2,_d2to40[i]);
    }
  }
  PORTABLE_INLINE_FUNCTION void
  DensityEnergyFromPressureTemperature(const Real press, const Real temp, Real *lambda,
                                       Real &rho, Real &sie) const;
  inline void Finalize() {}
  static std::string EosType() { return std::string("Vinet"); }
  static std::string EosPyType() { return EosType(); }

 private:
  static constexpr const unsigned long _preferred_input =
      thermalqs::density | thermalqs::specific_internal_energy;
  Real _rho0, _T0, _B0, _BP0, _A0, _Cv0, _E0, _S0;
  const Real *_expconsts;
  // static constexpr const char _eos_type[] = {"Vinet"};
  Real _VIP[42], _d2to40[38];

  PORTABLE_INLINE_FUNCTION void CheckVinet();
  PORTABLE_INLINE_FUNCTION void Vinet_F_DT_func(const Real rho, const Real T, 
		  				Real *output) const;
};

PORTABLE_INLINE_FUNCTION void Vinet::CheckVinet() {

  int ind;

  if (_rho0 < 0.0){
    printf("testing rho0\n");
    EOS_ERROR("Required Vinet model parameter rho0 not set");
    return;
  }
  if (_T0 < 0.0){
    printf("testing T0\n");
    EOS_ERROR("Required Vinet model parameter T0 not set");
    return;
  }
  if (_B0 < 0.0){
    printf("testing B0\n");
    EOS_ERROR("Required Vinet model parameter B0 not set");
    return;
  }
  if (_BP0 < 1.0){
    printf("testing BP0\n");
    EOS_ERROR("Required Vinet model parameter BP0 not set or too small (<1)");
    return;
  }
  if (_A0 < 0.0){
    printf("testing A0\n");
    EOS_ERROR("Required Vinet model parameter A0 not set");
    return;
  }
  if (_Cv0 < 0.0){
    printf("testing Cv0\n");
    EOS_ERROR("Required Vinet model parameter Cv0 not set");
    return;
  }


  for(ind = 0; ind < 38; ind++) _d2to40[ind]=_expconsts[ind];

  _VIP[0] = _B0/_rho0 + pow(_A0*_B0/_rho0,2.0)*_T0/_Cv0;      // sound speed squared
  _VIP[1] = 3.0/2.0*(_BP0-1.0); // exponent
  _VIP[2] = 1.0; // prefactor f0
  _VIP[3] = 0.0; // prefactor f1

  for (ind = 4; ind < 42; ind++) _VIP[ind] = _d2to40[ind-4];
  ind = 42;
  while (_VIP[ind] == 0.0) ind--;
  while (ind >= 4) {
    ind--;
    _VIP[ind] = _VIP[ind] - (ind -1)/_VIP[1]*_VIP[ind+1];
  }
}


PORTABLE_INLINE_FUNCTION void Vinet::Vinet_F_DT_func(const Real rho, const Real T, 
							Real *output) const {
  int pref0vp = -2, pref0vip = 2, ind;
  Real x, x2inv, onemx, etatimes1mx, expetatimes1mx, temp;
  Real sumP=0.0, sumB=0.0, sumE=0.0;
  Real entropy, energy, pressure, dpdrho, dpdt, dedt, dedrho, soundspeed;

  ind = 40-1;
  while (_VIP[pref0vip+ind] == 0.0) ind--;

  x = pow(_rho0/rho,1.0/3.0);  /*rho-dependent*/
  x2inv = 1.0/pow(x,2.0);
  onemx = 1.0-x;
  etatimes1mx = _VIP[1]*onemx;
  expetatimes1mx = exp(etatimes1mx);

  while (ind >= 2){
    sumP = _d2to40[pref0vp+ind]+onemx*sumP;
    sumB = _d2to40[pref0vp+ind]*ind+onemx*sumB;
    sumE = _VIP[pref0vip+ind]+onemx*sumE;
    ind--;
  }
  while (ind >= 0){
    sumE = _VIP[pref0vip+ind]+onemx*sumE;
    ind--;
  }
  sumP = 1.0+sumP*pow(onemx,2.0);
  sumB = sumB*onemx;
  sumE = _VIP[1]*onemx*sumE;

  /* T0 isotherm */
  energy = 9.0*_B0/pow(_VIP[1],2.0)/_rho0*(_VIP[pref0vip]-expetatimes1mx*(_VIP[pref0vip]-sumE))-_A0*_B0*_T0*(1.0-_rho0/rho);
  pressure = 3.0*_B0*x2inv*onemx*expetatimes1mx*sumP;
  temp = (1.0+onemx*(_VIP[1]*x+1.0))*sumP+x*onemx*sumB;
  temp = _B0/_rho0*x*expetatimes1mx*temp;

  /* Go to required temperature */
  energy = energy + _Cv0*(T-_T0)+_E0;
  pressure = pressure + _A0*_B0*(T-_T0);
  dpdrho = temp;
  dpdt = _A0*_B0;
  dedt = _Cv0;
  dedrho = (pressure - T*_A0*_B0)/pow(rho,2.0);
  if (T < 0.0){
    printf("warning!! negative temperature");
    entropy = 0.0;
  } else
    entropy = _A0*_B0*(1.0/rho-1/_rho0)+_Cv0*log(T/_T0)+_S0;
  soundspeed = std::max(0.0,dpdrho+T*pow(dpdt/rho,2.0)/_Cv0);
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
PORTABLE_INLINE_FUNCTION Real Vinet::InternalEnergyFromDensityTemperature(
    const Real rho, const Real temp, Real *lambda) const {
  Real output[8];
  Vinet_F_DT_func(rho,temp,output);
  return output[0];
}
PORTABLE_INLINE_FUNCTION Real Vinet::PressureFromDensityTemperature(
    const Real rho, const Real temp, Real *lambda) const {
  Real output[8];
  Vinet_F_DT_func(rho,temp,output);
  return output[1];
}
PORTABLE_INLINE_FUNCTION Real Vinet::EntropyFromDensityTemperature(
    const Real rho, const Real temp, Real *lambda) const {
  Real output[8];
  Vinet_F_DT_func(rho,temp,output);
  return output[6];
}
PORTABLE_INLINE_FUNCTION Real Vinet::BulkModulusFromDensityTemperature(
    const Real rho, const Real temp, Real *lambda) const {
  Real output[8];
  Vinet_F_DT_func(rho,temp,output);
  return output[7]*output[7]*rho;
}
PORTABLE_INLINE_FUNCTION Real Vinet::TemperatureFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  Real Tref;
  Real output[8];
  Tref=_T0;
  Vinet_F_DT_func(rho,Tref,output);
  return (sie - output[0])/_Cv0+Tref;
}
PORTABLE_INLINE_FUNCTION Real Vinet::PressureFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  Real temp;
  Real output[8];
  temp=TemperatureFromDensityInternalEnergy(rho,sie);
  Vinet_F_DT_func(rho,temp,output);
  return output[1];
}
PORTABLE_INLINE_FUNCTION Real Vinet::EntropyFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  Real temp;
  Real output[8];
  temp=TemperatureFromDensityInternalEnergy(rho,sie);
  Vinet_F_DT_func(rho,temp,output);
  return output[6];
}
PORTABLE_INLINE_FUNCTION Real Vinet::BulkModulusFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  Real temp;
  Real output[8];
  temp=TemperatureFromDensityInternalEnergy(rho,sie);
  Vinet_F_DT_func(rho,temp,output);
  return output[7]*output[7]*rho;
}
// Below are "unimplemented" routines
PORTABLE_INLINE_FUNCTION void Vinet::DensityEnergyFromPressureTemperature(
    const Real press, const Real temp, Real *lambda, Real &rho, Real &sie) const {
      EOS_ERROR("Vinet::DensityEnergyFromPressureTemperature: "
                "Not implemented.\n");
}
PORTABLE_INLINE_FUNCTION void Vinet::FillEos(Real &rho, Real &temp, Real &sie,
                                                 Real &press, Real &cv, Real &bmod,
                                                 const unsigned long output,
                                                 Real *lambda) const {
  // The following could be sped up with work!
  const unsigned long input = ~output;
  if (thermalqs::temperature & input && thermalqs::pressure & input) {
    DensityEnergyFromPressureTemperature(press, temp, lambda, rho, sie);
  } else if (thermalqs::density & output ||
             thermalqs::specific_internal_energy & output) {
    EOS_ERROR("Vinet FillEos: Density and energy are currently required inputs "
              "except when pressure and temperature are inputs.\n");
  }
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
  Real entropy;

  rho = _rho0;
  temp = _T0;
  sie = _E0;
  press = PressureFromDensityTemperature(rho, temp, lambda);
  entropy = _S0;
  cv = _Cv0;
  bmod = _B0;
  dpde = _A0*bmod/_Cv0;
  dvdt = _A0/_rho0;
}
} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_VINET_HPP_

//------------------------------------------------------------------------------
// © 2021-2023. Triad National Security, LLC. All rights reserved.  This
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

#include <ports-of-call/portable_errors.hpp>
#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/math_utils.hpp>
#include <singularity-eos/base/robust_utils.hpp>
#include <singularity-eos/eos/eos_base.hpp>

namespace singularity {

using namespace eos_base;
using namespace robust;

// COMMENT: This is a complete, thermodynamically consistent Vinet
// EOS created at Sandia National Laboratories by AEM.
// J. Appl. Phys. 119, 015904 (2016)
// It is particularly suited for high-pressure applications.
class Vinet : public EosBase<Vinet> {
 public:
  Vinet() = default;
  // Constructor
  Vinet(const Real rho0, const Real T0, const Real B0, const Real BP0, const Real A0,
        const Real Cv0, const Real E0, const Real S0, const Real *expconsts)
      : _rho0(rho0), _T0(T0), _B0(B0), _BP0(BP0), _A0(A0), _Cv0(Cv0), _E0(E0), _S0(S0) {
    CheckVinet();
    InitializeVinet(expconsts);
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
  PORTABLE_INLINE_FUNCTION Real
  MinInternalEnergyFromDensity(const Real rho, Real *lambda = nullptr) const;
  // Entropy added AEM Dec. 2022
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
  // Thermal Bulk Modulus added AEM Dec 2022
  PORTABLE_INLINE_FUNCTION Real TBulkModulusFromDensityTemperature(
      const Real rho, const Real temp, Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityTemperature(
      const Real rho, const Real temp, Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  // Thermal expansion coefficient added AEM 2022
  PORTABLE_INLINE_FUNCTION Real TExpansionCoeffFromDensityTemperature(
      const Real rho, const Real temp, Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temp, Real *lambda = nullptr) const {
    return robust::ratio(_A0 * _B0, _Cv0 * rho);
  }
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const {
    return robust::ratio(_A0 * _B0, _Cv0 * rho);
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
  static inline unsigned long scratch_size(std::string method, unsigned int nelements) {
    return 0;
  }
  static inline unsigned long max_scratch_size(unsigned int nelements) { return 0; }
  PORTABLE_INLINE_FUNCTION void PrintParams() const {
    static constexpr char st[]{"Vinet Params: "};
    printf("%s rho0:%e T0:%e B0:%e BP0:%e\n  A0:%e Cv0:%e E0:%e S0:%e\n"
           "non-zero elements in d2tod40 array:\n",
           st, _rho0, _T0, _B0, _BP0, _A0, _Cv0, _E0, _S0);
    for (int i = 0; i < 39; i++) {
      if (_d2tod40[i] > 0.0) printf("d%i:%e\t", i + 2, _d2tod40[i]);
    }
    printf("\n\n");
  }
  // Density/Energy from P/T not unique, if used will give error
  PORTABLE_INLINE_FUNCTION void
  DensityEnergyFromPressureTemperature(const Real press, const Real temp, Real *lambda,
                                       Real &rho, Real &sie) const;
  inline void Finalize() {}
  static std::string EosType() { return std::string("Vinet"); }
  static std::string EosPyType() { return EosType(); }

 private:
  static constexpr const unsigned long _preferred_input =
      thermalqs::density | thermalqs::temperature;
  Real _rho0, _T0, _B0, _BP0, _A0, _Cv0, _E0, _S0;
  static constexpr const int PressureCoeffsd2tod40Size = 39;
  static constexpr const int VinetInternalParametersSize = PressureCoeffsd2tod40Size + 4;
  Real _VIP[VinetInternalParametersSize], _d2tod40[PressureCoeffsd2tod40Size];
  void CheckVinet();
  void InitializeVinet(const Real *expcoeffs);
  PORTABLE_INLINE_FUNCTION void Vinet_F_DT_func(const Real rho, const Real T,
                                                Real *output) const;
};

inline void Vinet::CheckVinet() {

  if (_rho0 < 0.0) {
    PORTABLE_ALWAYS_THROW_OR_ABORT("Required Vinet model parameter rho0 < 0");
  }
  if (_T0 < 0.0) {
    PORTABLE_ALWAYS_THROW_OR_ABORT("Required Vinet model parameter T0 < 0");
  }
  if (_B0 < 0.0) {
    PORTABLE_ALWAYS_THROW_OR_ABORT("Required Vinet model parameter B0 < 0");
  }
  if (_BP0 < 1.0) {
    PORTABLE_ALWAYS_THROW_OR_ABORT("Required Vinet model parameter BP0 < 1");
  }
  if (_A0 < 0.0) {
    PORTABLE_ALWAYS_THROW_OR_ABORT("Required Vinet model parameter A0 < 0");
  }
  if (_Cv0 < 0.0) {
    PORTABLE_ALWAYS_THROW_OR_ABORT("Required Vinet model parameter Cv0 < 0");
  }
}

inline void Vinet::InitializeVinet(const Real *d2tod40input) {

  // The PressureCoeffsd2tod40Size (=39) allowed d2 to d40 coefficients
  // for the pressure reference curve vs rho
  // are seldom all used so did not want to crowd the argument list with them.
  // Instead I ask the host code to send me a pointer to this array so that I can
  // copy it here. Not used coeffs should be set to 0.0 (of course).
  for (int ind = 0; ind < PressureCoeffsd2tod40Size; ind++) {
    _d2tod40[ind] = d2tod40input[ind];
  }
  // Put a couple of  much used model parameter combinations and
  // the energy coefficients f0 to f40 in an internal parameters array
  // of size 2+41=VinetInternalParametersSize.
  _VIP[0] = robust::ratio(_B0, _rho0) +
            robust::ratio(_A0 * _A0 * _B0 * _B0 * _T0,
                          _rho0 * _rho0 * _Cv0); // sound speed squared
  _VIP[1] = 3.0 / 2.0 * (_BP0 - 1.0);            // exponent eta0
  // initializing with pressure coeffs to get energy coeffs into VIP
  _VIP[2] = 1.0;                                                // prefactor d0
  _VIP[3] = 0.0;                                                // prefactor d1
  for (int ind = 4; ind < VinetInternalParametersSize; ind++) { //_d2tod40[0]=d2
    _VIP[ind] = _d2tod40[ind - 4];                              // dn is in VIP[n+2]
  }
  for (int ind = VinetInternalParametersSize - 2; ind >= 2;
       ind--) { // _VIP[42]=d40=f40 given,first calculated is _VIP[41]=f39
    _VIP[ind] = _VIP[ind] - (ind) / _VIP[1] * _VIP[ind + 1]; // prefactors f40 to f0
  }                                                          // _VIP[n+2]=fn, ind=n+2
}

PORTABLE_INLINE_FUNCTION void Vinet::Vinet_F_DT_func(const Real rho, const Real T,
                                                     Real *output) const {
  constexpr int pref0vp = -2, pref0vip = 2, maxind = VinetInternalParametersSize - 3;
  Real sumP = 0.0, sumB = 0.0, sumE = 0.0;

  Real x = std::cbrt(robust::ratio(_rho0, rho)); /*rho-dependent*/
  Real x2inv = robust::ratio(1.0, x * x);
  Real onemx = 1.0 - x;
  Real etatimes1mx = _VIP[1] * onemx;
  Real expetatimes1mx = exp(etatimes1mx);

#pragma unroll
  for (int ind = maxind; ind >= 2; ind--) {              //_d2tod40[0]=d2
    sumP = _d2tod40[pref0vp + ind] + onemx * sumP;       //_d2tod40[38]=d40
    sumB = _d2tod40[pref0vp + ind] * ind + onemx * sumB; //_d2tod40[-2+40]=d40
    sumE = _VIP[pref0vip + ind] + onemx * sumE;          //_VIP[42]=f40
  }                                                      //_VIP[2]=f0
#pragma unroll
  for (int ind = 1; ind >= 0; ind--) {
    sumE = _VIP[pref0vip + ind] + onemx * sumE;
  }
  sumP = 1.0 + sumP * (onemx * onemx);
  sumB = sumB * onemx;
  sumE = _VIP[1] * onemx * sumE;

  /* T0 isotherm */
  Real energy = 9.0 * robust::ratio(_B0, _VIP[1] * _VIP[1] * _rho0) *
                    (_VIP[pref0vip] - expetatimes1mx * (_VIP[pref0vip] - sumE)) -
                (_A0 * _B0) * (robust::ratio(_T0, _rho0) - robust::ratio(_T0, rho));
  Real pressure = 3.0 * _B0 * x2inv * onemx * expetatimes1mx * sumP;
  Real temp = (1.0 + onemx * (_VIP[1] * x + 1.0)) * sumP + x * onemx * sumB;
  temp = robust::ratio(_B0, _rho0) * x * expetatimes1mx * temp;

  /* Go to required temperature */
  energy = energy + _Cv0 * (T - _T0) + _E0;
  pressure = pressure + (_A0 * _B0) * (T - _T0);
  Real dpdrho = temp;
  Real dpdt = _A0 * _B0;
  Real dedt = _Cv0;
  Real dedrho = robust::ratio(pressure - T * (_A0 * _B0), rho * rho);
  Real entropy;
  if (T < 0.0) {
#ifndef NDEBUG
    PORTABLE_WARN("Negative temperature input");
#endif // NDEBUG
    entropy = 0.0;
  } else {
    entropy = (_A0 * _B0) * (robust::ratio(1.0, rho) - robust::ratio(1, _rho0)) +
              _Cv0 * std::log(robust::ratio(T, _T0)) + _S0;
  }
  Real soundspeed =
      std::max(0.0, dpdrho + T * robust::ratio(dpdt * dpdt, rho * rho * _Cv0));
  soundspeed = std::sqrt(soundspeed);

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
  Vinet_F_DT_func(rho, temp, output);
  return output[0];
}
PORTABLE_INLINE_FUNCTION Real Vinet::PressureFromDensityTemperature(const Real rho,
                                                                    const Real temp,
                                                                    Real *lambda) const {
  Real output[8];
  Vinet_F_DT_func(rho, temp, output);
  return output[1];
}
PORTABLE_INLINE_FUNCTION Real Vinet::EntropyFromDensityTemperature(const Real rho,
                                                                   const Real temp,
                                                                   Real *lambda) const {
  Real output[8];
  Vinet_F_DT_func(rho, temp, output);
  return output[6];
}
PORTABLE_INLINE_FUNCTION Real Vinet::TExpansionCoeffFromDensityTemperature(
    const Real rho, const Real temp, Real *lambda) const {
  Real output[8];
  Vinet_F_DT_func(rho, temp, output);
  return robust::ratio(output[3], output[2] * rho);
}
PORTABLE_INLINE_FUNCTION Real Vinet::TBulkModulusFromDensityTemperature(
    const Real rho, const Real temp, Real *lambda) const {
  Real output[8];
  Vinet_F_DT_func(rho, temp, output);
  return output[2] * rho;
}
PORTABLE_INLINE_FUNCTION Real Vinet::BulkModulusFromDensityTemperature(
    const Real rho, const Real temp, Real *lambda) const {
  Real output[8];
  Vinet_F_DT_func(rho, temp, output);
  return output[7] * output[7] * rho;
}
PORTABLE_INLINE_FUNCTION Real Vinet::TemperatureFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  Real Tref;
  Real output[8];
  Tref = _T0;
  Vinet_F_DT_func(rho, Tref, output);
  return robust::ratio(sie - output[0], _Cv0) + Tref;
}
PORTABLE_INLINE_FUNCTION Real Vinet::PressureFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  Real temp;
  Real output[8];
  temp = TemperatureFromDensityInternalEnergy(rho, sie);
  Vinet_F_DT_func(rho, temp, output);
  return output[1];
}
PORTABLE_INLINE_FUNCTION Real Vinet::MinInternalEnergyFromDensity(const Real rho,
                                                                  Real *lambda) const {
#ifndef NDEBUG
      printf(
          "WARNING: MinInternalEnergtyFromDensity is not defined for Vinet EOS.");
#endif

  return 0.0;
}
PORTABLE_INLINE_FUNCTION Real Vinet::EntropyFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  Real temp;
  Real output[8];
  temp = TemperatureFromDensityInternalEnergy(rho, sie);
  Vinet_F_DT_func(rho, temp, output);
  return output[6];
}
PORTABLE_INLINE_FUNCTION Real Vinet::BulkModulusFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  Real temp;
  Real output[8];
  temp = TemperatureFromDensityInternalEnergy(rho, sie);
  Vinet_F_DT_func(rho, temp, output);
  return output[7] * output[7] * rho;
}
// AEM: Give error since function is not well defined
PORTABLE_INLINE_FUNCTION void
Vinet::DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                            Real *lambda, Real &rho, Real &sie) const {
  EOS_ERROR("Vinet::DensityEnergyFromPressureTemperature: "
            "Not implemented.\n");
}
// AEM: We should add entropy and Gruneissen parameters here so that it is complete
// If we add also alpha and BT, those should also be in here.
PORTABLE_INLINE_FUNCTION void Vinet::FillEos(Real &rho, Real &temp, Real &sie,
                                             Real &press, Real &cv, Real &bmod,
                                             const unsigned long output,
                                             Real *lambda) const {
  const unsigned long input = ~output; // everything that is not output is input
  if (thermalqs::density & output) {
    EOS_ERROR("Vinet FillEos: Density is required input.\n");
  } else if (thermalqs::temperature & output &&
             thermalqs::specific_internal_energy & output) {
    EOS_ERROR("Vinet FillEos: Density and Internal Energy or Density and Temperature "
              "are required input parameters.\n");
  }
  if (thermalqs::specific_internal_energy & input) {
    temp = TemperatureFromDensityInternalEnergy(rho, sie);
  }
  Real Vout[8];
  Vinet_F_DT_func(rho, temp, Vout);
  if (output & thermalqs::temperature) temp = temp;
  if (output & thermalqs::specific_internal_energy) sie = Vout[0];
  if (output & thermalqs::pressure) press = Vout[1];
  if (output & thermalqs::specific_heat) cv = Vout[4];
  if (output & thermalqs::bulk_modulus) bmod = Vout[7] * Vout[7] * rho;
}

// TODO(JMM): pre-cache these rather than recomputing them each time
PORTABLE_INLINE_FUNCTION
void Vinet::ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press,
                                   Real &cv, Real &bmod, Real &dpde, Real &dvdt,
                                   Real *lambda) const {
  // AEM: Added all variables I think should be output eventually
  Real tbmod;
  // Real entropy, alpha, Gamma;

  rho = _rho0;
  temp = _T0;
  sie = _E0;
  press = PressureFromDensityTemperature(rho, temp, lambda);
  // entropy = _S0;
  cv = _Cv0;
  tbmod = _B0;
  // alpha = _A0;
  bmod = BulkModulusFromDensityTemperature(rho, temp, lambda);
  // Gamma = robust::ratio(_A0 * _B0, _Cv0 * _rho0);
  // AEM: I suggest taking the two following away.
  dpde = robust::ratio(_A0 * tbmod, _Cv0);
  dvdt = robust::ratio(_A0, _rho0);
}
} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_VINET_HPP_

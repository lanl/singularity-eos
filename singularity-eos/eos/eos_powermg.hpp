//------------------------------------------------------------------------------
// © 2021-2024. Triad National Security, LLC. All rights reserved.  This
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

#ifndef _SINGULARITY_EOS_EOS_POWERMG_HPP_
#define _SINGULARITY_EOS_EOS_POWERMG_HPP_

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

// COMMENT: This is a complete, thermodynamically consistent Mie-Gruneisen
// EOS based on a power series in eta for the Hugoniot pressure.
class PowerMG : public EosBase<PowerMG> {
 public:
  PowerMG() = default;
  // Constructor
  PowerMG(const Real rho0, const Real T0, const Real G0, const Real Cv0, const Real E0,
          const Real S0, const Real Pmin, const Real *expconsts)
      : _rho0(rho0), _T0(T0), _G0(G0), _Cv0(Cv0), _E0(E0), _S0(S0), _Pmin(Pmin) {
    _InitializePowerMG(expconsts);
    _CheckPowerMG();
  }

  PowerMG GetOnDevice() { return *this; }
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
  // added for testing AEM Dec 2023
  PORTABLE_INLINE_FUNCTION Real AllHugPressureFromDensity(const Real rho) const;
  PORTABLE_INLINE_FUNCTION Real AllHugInternalEnergyFromDensity(const Real rho) const;
  PORTABLE_INLINE_FUNCTION Real AllHugTemperatureFromDensity(const Real rho) const;
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
    return robust::ratio(_G0 * _rho0, rho);
  }
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const {
    return robust::ratio(_G0 * _rho0, rho);
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
  SG_ADD_BASE_CLASS_USINGS(PowerMG)
  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return 0; }
  static constexpr unsigned long PreferredInput() { return _preferred_input; }
  static inline unsigned long scratch_size(std::string method, unsigned int nelements) {
    return 0;
  }
  static inline unsigned long max_scratch_size(unsigned int nelements) { return 0; }
  PORTABLE_INLINE_FUNCTION void PrintParams() const {
    static constexpr char st[]{"PowerMG Params: "};
    printf("%s rho0:%e T0:%e G0:%e Cv0:%e E0:%e S0:%e Pmin:%e\n"
           "non-zero elements in K0toK40 array:\n",
           st, _rho0, _T0, _G0, _Cv0, _E0, _S0, _Pmin);
    for (int i = 0; i < 41; i++) {
      if (_K0toK40[i] * _K0toK40[i] > 0.0) printf("K%i:%e\t", i, _K0toK40[i]);
    }
    printf("\n\n");
  }
  // Density/Energy from P/T not unique, if used will give error
  PORTABLE_INLINE_FUNCTION void
  DensityEnergyFromPressureTemperature(const Real press, const Real temp, Real *lambda,
                                       Real &rho, Real &sie) const;
  inline void Finalize() {}
  static std::string EosType() { return std::string("PowerMG"); }
  static std::string EosPyType() { return EosType(); }

 private:
  static constexpr const unsigned long _preferred_input =
      thermalqs::density | thermalqs::specific_internal_energy;
  Real _rho0, _T0, _G0, _Cv0, _E0, _S0, _Pmin;
  static constexpr const int PressureCoeffsK0toK40Size = 41;
  int _M;
  Real _K0toK40[PressureCoeffsK0toK40Size];
  void _CheckPowerMG();
  void _InitializePowerMG(const Real *expcoeffs);
  PORTABLE_INLINE_FUNCTION Real _HugPressureFromDensity(const Real rho) const;
  PORTABLE_INLINE_FUNCTION Real _HugTemperatureFromDensity(const Real rho) const;
  PORTABLE_INLINE_FUNCTION Real _SN2Mp2(const Real x) const;
  PORTABLE_INLINE_FUNCTION Real
  _compBulkModulusFromDensityTemperature(const Real rho, const Real temp) const;
  PORTABLE_INLINE_FUNCTION Real
  _compBulkModulusFromDensityInternalEnergy(const Real rho, const Real sie) const;
};

inline void PowerMG::_CheckPowerMG() {

  if (_rho0 < 0.0) {
    PORTABLE_ALWAYS_THROW_OR_ABORT("Required PowerMG model parameter rho0 < 0");
  }
  if (_T0 < 0.0) {
    PORTABLE_ALWAYS_THROW_OR_ABORT("Required PowerMG model parameter T0 < 0");
  }
  if (_Cv0 < 0.0) {
    PORTABLE_ALWAYS_THROW_OR_ABORT("Required PowerMG model parameter Cv0 < 0");
  }
  if (_K0toK40[0] < 0.0) {
    PORTABLE_ALWAYS_THROW_OR_ABORT("Required PowerMG model parameter K0 < 0");
  }
  if (_Pmin >= 0.0) {
    _Pmin = -1000 * _K0toK40[0];
#ifndef NDEBUG
    PORTABLE_WARN(
        "PowerMG model parameter Pmin not set or positive. Reset to default (-1000*K0)");
#endif // NDEBUG
  }
}

inline void PowerMG::_InitializePowerMG(const Real *K0toK40input) {

  // The PressureCoeffsK0toK40Size (=41) allowed K0 to K40 coefficients
  // for the pressure reference curve vs rho
  // are seldom all used so did not want to crowd the argument list with them.
  // Instead I ask the host code to send me a pointer to this array so that I can
  // copy it here. Not used coeffs should be set to 0.0 (of course).
  _M = 0;
  for (int ind = 0; ind < PressureCoeffsK0toK40Size; ind++) {
    _K0toK40[ind] = K0toK40input[ind];
    if (_K0toK40[ind] != 0.0) _M = ind;
  }
}

PORTABLE_INLINE_FUNCTION Real PowerMG::_HugPressureFromDensity(Real rho) const {
  const Real eta = 1.0 - robust::ratio(_rho0, rho);
  Real value = 0.0;
  for (int ind = _M; ind >= 1; ind--) {
    value = eta * value + _K0toK40[ind];
  }
  value = _K0toK40[0] * eta * (eta * value + 1.0);
  return value;
}
PORTABLE_INLINE_FUNCTION Real PowerMG::_HugTemperatureFromDensity(Real rho) const {
  const Real eta = 1.0 - robust::ratio(_rho0, rho);
  Real gamma2 = _SN2Mp2(_G0 * eta);
  Real sum = _M * _K0toK40[_M] * gamma2;
  for (int ind = _M - 1; ind >= 1; ind--) {
    gamma2 = _G0 * eta / (ind + 2) * gamma2 + 1.0 / (ind + 2);
    sum = eta * sum + ind * _K0toK40[ind] * gamma2;
  }
  Real temp = eta * eta * eta * sum * (_K0toK40[0] / _Cv0 / 2.0 / _rho0);
  temp = _T0 * exp(eta * _G0) + temp;
  return temp;
}
// The function S_2^N is described in Robinson's report (SAND2019-6025
// equations 1.29-3nd 1.30 with n=M+2)  and is the largest coefficient in
// the series in eta for the Hugoniot temperature, equation 1.34. Other
// coefficients are calculated recusively from this one according to equation 1.36.
PORTABLE_INLINE_FUNCTION Real PowerMG::_SN2Mp2(const Real x) const {
  Real ind = _M + 2;
  int maxind = 200;
  Real ak = 1.0 / ind;
  Real sum = ak;
  Real pf = exp(x) - 1.0;
  Real ek = pf * ak;
  // Add terms to the series sum until the relative error ek
  // is smaller than machine presision.
  while ((ek > sum * 1.e-15) && (ind < maxind)) {
    ind = ind + 1;
    ak = ak * x / ind;
    ek = pf * ak;
    sum = sum + ak;
  }
  if (ind >= maxind) {
#ifndef NDEBUG
    PORTABLE_WARN("Part of Hugoniot Temperature not converged");
#endif // NDEBUG
  }
  return sum;
}

PORTABLE_INLINE_FUNCTION Real PowerMG::AllHugPressureFromDensity(Real rho) const {
  const Real eta = 1.0 - robust::ratio(_rho0, rho);
  Real value = 0.0;
  if (eta <= 0.0) {
    if (eta < _Pmin / _K0toK40[0]) {
      value = _Pmin;
    } else {
      value = _K0toK40[0] * eta;
    }
  } else {
    value = _HugPressureFromDensity(rho);
  }
  return value;
}
PORTABLE_INLINE_FUNCTION Real PowerMG::AllHugInternalEnergyFromDensity(Real rho) const {
  const Real eta = 1.0 - robust::ratio(_rho0, rho);
  return _E0 + robust::ratio(eta, 2.0 * _rho0) * AllHugPressureFromDensity(rho);
}
PORTABLE_INLINE_FUNCTION Real PowerMG::AllHugTemperatureFromDensity(Real rho) const {
  const Real eta = 1.0 - robust::ratio(_rho0, rho);
  Real value = 0.0;
  if (eta <= 0.0) {
    value = _T0 * exp(eta * _G0);
  } else {
    value = _HugTemperatureFromDensity(rho);
  }
  return value;
}

PORTABLE_INLINE_FUNCTION Real PowerMG::InternalEnergyFromDensityTemperature(
    const Real rho, const Real temp, Real *lambda) const {
  Real value = AllHugInternalEnergyFromDensity(rho) +
               _Cv0 * (temp - AllHugTemperatureFromDensity(rho));
  return value;
}
PORTABLE_INLINE_FUNCTION Real PowerMG::PressureFromDensityTemperature(
    const Real rho, const Real temp, Real *lambda) const {
  Real value = AllHugPressureFromDensity(rho) +
               _G0 * _rho0 * _Cv0 * (temp - AllHugTemperatureFromDensity(rho));
  return value;
}
PORTABLE_INLINE_FUNCTION Real PowerMG::EntropyFromDensityTemperature(const Real rho,
                                                                     const Real temp,
                                                                     Real *lambda) const {
  const Real eta = 1.0 - robust::ratio(_rho0, rho);
  Real value = _S0 - _G0 * _Cv0 * eta + _Cv0 * std::log(temp / _T0);
  return value;
}
PORTABLE_INLINE_FUNCTION Real PowerMG::TExpansionCoeffFromDensityTemperature(
    const Real rho, const Real temp, Real *lambda) const {
  Real value =
      robust::ratio(_Cv0 * _rho0 * _G0, TBulkModulusFromDensityTemperature(rho, temp));
  return value;
}
PORTABLE_INLINE_FUNCTION Real PowerMG::TBulkModulusFromDensityTemperature(
    const Real rho, const Real temp, Real *lambda) const {
  return BulkModulusFromDensityTemperature(rho, temp) -
         _G0 * _G0 * _Cv0 * _rho0 * robust::ratio(_rho0, rho) * temp;
}
PORTABLE_INLINE_FUNCTION Real
PowerMG::_compBulkModulusFromDensityTemperature(const Real rho, const Real temp) const {
  const Real eta = 1.0 - robust::ratio(_rho0, rho);
  Real value = (1.0 - _G0 * eta / 2.0) * _K0toK40[0];
  Real sum = 0.0;
  for (int ind = _M; ind >= 1; ind--) {
    sum = sum * eta + ind * _K0toK40[ind] * eta;
  }
  value = value * sum;
  value = value + robust::ratio(_HugPressureFromDensity(rho), eta);
  value = value + _G0 * _G0 * _Cv0 * _rho0 * (temp - _HugTemperatureFromDensity(rho));
  value = robust::ratio(_rho0, rho) * value;
  return value;
}
PORTABLE_INLINE_FUNCTION Real PowerMG::BulkModulusFromDensityTemperature(
    const Real rho, const Real temp, Real *lambda) const {
  const Real eta = 1.0 - robust::ratio(_rho0, rho);
  Real value;
  if (eta <= 0.0) {
    if (eta < _Pmin / _K0toK40[0]) {
      value = _G0 * _G0 * _Cv0 * _rho0 * (temp - _T0 * exp(_G0 * eta));
    } else {
      value = _K0toK40[0] + _G0 * _G0 * _Cv0 * _rho0 * (temp - _T0 * exp(_G0 * eta));
    }
    value = robust::ratio(_rho0, rho) * value;
  } else {
    value = _compBulkModulusFromDensityTemperature(rho, temp);
  }
  return value;
}
PORTABLE_INLINE_FUNCTION Real
PowerMG::_compBulkModulusFromDensityInternalEnergy(const Real rho, const Real sie) const {
  const Real eta = 1.0 - robust::ratio(_rho0, rho);
  Real value = (1.0 - _G0 * eta / 2.0) * _K0toK40[0];
  Real sum = 0.0;
  for (int ind = _M; ind >= 1; ind--) {
    sum = sum * eta + ind * _K0toK40[ind] * eta;
  }
  value = value * sum;
  value = value + robust::ratio(_HugPressureFromDensity(rho), eta);
  value = value + _G0 * _G0 * _rho0 * (sie - AllHugInternalEnergyFromDensity(rho));
  value = robust::ratio(_rho0, rho) * value;
  return value;
}
PORTABLE_INLINE_FUNCTION Real PowerMG::BulkModulusFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  const Real eta = 1.0 - robust::ratio(_rho0, rho);
  Real value;
  if (eta <= 0.0) {
    if (eta < _Pmin / _K0toK40[0]) {
      value = _G0 * _G0 * _rho0 *
              (sie - _E0 + _Pmin / _rho0 * (_Pmin / _K0toK40[0] / 2.0 - eta));
    } else {
      value =
          _K0toK40[0] + _G0 * _G0 * _rho0 * (sie - _E0 - _K0toK40[0] * eta * eta / 2.0);
    }
    value = robust::ratio(_rho0, rho) * value;
  } else {
    value = _compBulkModulusFromDensityInternalEnergy(rho, sie);
  }
  return value;
}

PORTABLE_INLINE_FUNCTION Real PowerMG::TemperatureFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  Real value = (sie - AllHugInternalEnergyFromDensity(rho)) / _Cv0 +
               AllHugTemperatureFromDensity(rho);
  if (value < 0.0) {
#ifndef NDEBUG
    PORTABLE_WARN("Negative temperature");
#endif // NDEBUG
    value = 1.e-12;
  }
  return value;
}
PORTABLE_INLINE_FUNCTION Real PowerMG::PressureFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  Real value = AllHugPressureFromDensity(rho) +
               _rho0 * _G0 * (sie - AllHugInternalEnergyFromDensity(rho));
  return value;
}
PORTABLE_INLINE_FUNCTION Real PowerMG::MinInternalEnergyFromDensity(const Real rho,
                                                                    Real *lambda) const {
  MinInternalEnergyIsNotEnabled("PowerMG");
  return 0.0;
}
PORTABLE_INLINE_FUNCTION Real PowerMG::EntropyFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  const Real eta = 1.0 - robust::ratio(_rho0, rho);
  Real value = std::log(TemperatureFromDensityInternalEnergy(rho, sie) / _T0);
  value = _S0 - _G0 * _Cv0 * eta + _Cv0 * value;
  if (value < 0.0) {
#ifndef NDEBUG
    PORTABLE_WARN("Negative entropy");
#endif // NDEBUG
    value = 1.e-12;
  }
  return value;
}
// AEM: Give error since function is not well defined
PORTABLE_INLINE_FUNCTION void
PowerMG::DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                              Real *lambda, Real &rho, Real &sie) const {
  EOS_ERROR("PowerMG::DensityEnergyFromPressureTemperature: "
            "Not implemented.\n");
}
// AEM: We should add entropy and Gruneissen parameters here so that it is complete
// If we add also alpha and BT, those should also be in here.
PORTABLE_INLINE_FUNCTION void PowerMG::FillEos(Real &rho, Real &temp, Real &sie,
                                               Real &press, Real &cv, Real &bmod,
                                               const unsigned long output,
                                               Real *lambda) const {
  const unsigned long input = ~output; // everything that is not output is input
  if (thermalqs::density & output) {
    EOS_ERROR("PowerMG FillEos: Density is required input.\n");
  } else if (thermalqs::temperature & output &&
             thermalqs::specific_internal_energy & output) {
    EOS_ERROR("PowerMG FillEos: Density and Internal Energy or Density and Temperature "
              "are required input parameters.\n");
  }
  if (thermalqs::temperature & input) {
    sie = InternalEnergyFromDensityTemperature(rho, temp);
  }
  if (output & thermalqs::temperature)
    temp = TemperatureFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::specific_internal_energy) sie = sie;
  if (output & thermalqs::pressure) press = PressureFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::specific_heat)
    cv = SpecificHeatFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::bulk_modulus)
    bmod = BulkModulusFromDensityInternalEnergy(rho, sie);
}

// TODO(JMM): pre-cache these rather than recomputing them each time
PORTABLE_INLINE_FUNCTION
void PowerMG::ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press,
                                     Real &cv, Real &bmod, Real &dpde, Real &dvdt,
                                     Real *lambda) const {
  // AEM: Added all variables I think should be output eventually
  //Real tbmod{};
  // Real entropy, alpha, Gamma;

  rho = _rho0;
  temp = _T0;
  sie = _E0;
  press = PressureFromDensityTemperature(rho, temp, lambda);
  // entropy = _S0;
  cv = _Cv0;
  //tbmod = _K0toK40[0] - _G0 * _G0 * _Cv0 * _rho0 * _T0;
  // alpha = _A0;
  bmod = _K0toK40[0];
  // Gamma = _G0;
  // AEM: I suggest taking the two following away.
  dpde = _G0 * _rho0;
  dvdt = robust::ratio(-1.0, _T0 * _G0 * _rho0);
}
} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_POWERMG_HPP_

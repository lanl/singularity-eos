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

#ifndef _SINGULARITY_EOS_EOS_MGUSUP_HPP_
#define _SINGULARITY_EOS_EOS_MGUSUP_HPP_

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
// EOS based on a linear Us-up Hugoniot, Us=Cs+s*up.
class MGUsup : public EosBase<MGUsup> {
 public:
  MGUsup() = default;
  // Constructor
  MGUsup(const Real rho0, const Real T0, const Real Cs, const Real s, const Real G0,
         const Real Cv0, const Real E0, const Real S0)
      : _rho0(rho0), _T0(T0), _Cs(Cs), _s(s), _G0(G0), _Cv0(Cv0), _E0(E0), _S0(S0) {
    _CheckMGUsup();
  }

  MGUsup GetOnDevice() { return *this; }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temp, Indexer_t &&lambda = nullptr) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityTemperature(
      const Real rho, const Real temp, Indexer_t &&lambda = nullptr) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real
  MinInternalEnergyFromDensity(const Real rho, Indexer_t &&lambda = nullptr) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real EntropyFromDensityTemperature(
      const Real rho, const Real temp, Indexer_t &&lambda = nullptr) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real EntropyFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityTemperature(
      const Real rho, const Real temp, Indexer_t &&lambda = nullptr) const {
    return _Cv0;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    return _Cv0;
  }
  // added for testing AEM Dec 2023
  PORTABLE_INLINE_FUNCTION Real HugPressureFromDensity(const Real rho) const;
  PORTABLE_INLINE_FUNCTION Real HugInternalEnergyFromDensity(const Real rho) const;
  PORTABLE_INLINE_FUNCTION Real HugTemperatureFromDensity(const Real rho) const;
  // Thermal Bulk Modulus added AEM Dec 2022
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real TBulkModulusFromDensityTemperature(
      const Real rho, const Real temp, Indexer_t &&lambda = nullptr) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityTemperature(
      const Real rho, const Real temp, Indexer_t &&lambda = nullptr) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const;
  // Thermal expansion coefficient added AEM 2022
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real TExpansionCoeffFromDensityTemperature(
      const Real rho, const Real temp, Indexer_t &&lambda = nullptr) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temp, Indexer_t &&lambda = nullptr) const {
    return robust::ratio(_G0 * _rho0, rho);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    return robust::ratio(_G0 * _rho0, rho);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void FillEos(Real &rho, Real &temp, Real &energy, Real &press,
                                        Real &cv, Real &bmod, const unsigned long output,
                                        Indexer_t &&lambda = nullptr) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                         Real &bmod, Real &dpde, Real &dvdt,
                         Indexer_t &&lambda = nullptr) const;
  // Generic functions provided by the base class. These contain e.g. the vector
  // overloads that use the scalar versions declared here
  SG_ADD_BASE_CLASS_USINGS(MGUsup)
  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return 0; }
  static constexpr unsigned long PreferredInput() { return _preferred_input; }
  static inline unsigned long scratch_size(std::string method, unsigned int nelements) {
    return 0;
  }
  static inline unsigned long max_scratch_size(unsigned int nelements) { return 0; }
  PORTABLE_INLINE_FUNCTION void PrintParams() const {
    static constexpr char st[]{"MGUsup Params: "};
    printf("%s rho0:%e T0:%e Cs:%e s:%e\n  G0:%e Cv0:%e E0:%e S0:%e\n", st, _rho0, _T0,
           _Cs, _s, _G0, _Cv0, _E0, _S0);
    printf("\n\n");
  }
  // Density/Energy from P/T not unique, if used will give error
  template <typename Indexer_t>
  PORTABLE_INLINE_FUNCTION void
  DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                       Indexer_t &&lambda, Real &rho, Real &sie) const;
  inline void Finalize() {}
  static std::string EosType() { return std::string("MGUsup"); }
  static std::string EosPyType() { return EosType(); }

 private:
  static constexpr const unsigned long _preferred_input =
      thermalqs::density | thermalqs::specific_internal_energy;
  Real _rho0, _T0, _Cs, _s, _G0, _Cv0, _E0, _S0;
  void _CheckMGUsup();
};

inline void MGUsup::_CheckMGUsup() {

  if (_rho0 < 0.0) {
    PORTABLE_ALWAYS_THROW_OR_ABORT("Required MGUsup model parameter rho0 < 0");
  }
  if (_T0 < 0.0) {
    PORTABLE_ALWAYS_THROW_OR_ABORT("Required MGUsup model parameter T0 < 0");
  }
  if (_Cs < 0.0) {
    PORTABLE_ALWAYS_THROW_OR_ABORT("Required MGUsup model parameter Cs < 0");
  }
  if (_s < 0.0) {
    PORTABLE_ALWAYS_THROW_OR_ABORT("Required MGUsup model parameter s < 0");
  }
  if (_G0 < 0.0) {
    PORTABLE_ALWAYS_THROW_OR_ABORT("Required MGUsup model parameter G0 < 0");
  }
  if (_Cv0 < 0.0) {
    PORTABLE_ALWAYS_THROW_OR_ABORT("Required MGUsup model parameter Cv0 < 0");
  }
}

PORTABLE_INLINE_FUNCTION Real MGUsup::HugPressureFromDensity(Real rho) const {
  Real eta = 1.0 - robust::ratio(_rho0, rho);
  return _Cs * _Cs * _rho0 * robust::ratio(eta, (1.0 - _s * eta) * (1.0 - _s * eta));
}
PORTABLE_INLINE_FUNCTION Real MGUsup::HugInternalEnergyFromDensity(Real rho) const {
  Real eta = 1.0 - robust::ratio(_rho0, rho);
  return _E0 + robust::ratio(eta, 2.0 * _rho0) * HugPressureFromDensity(rho);
}
PORTABLE_INLINE_FUNCTION Real MGUsup::HugTemperatureFromDensity(Real rho) const {
  int sumkmax = 20;
  Real cutoff = 1.e-16;
  Real eta = 1.0 - robust::ratio(_rho0, rho);
  Real f1 = 1.0 - _s * eta;
  if (f1 <= 0.0) {
    PORTABLE_ALWAYS_THROW_OR_ABORT("MGUsup model parameters s and rho0 together with rho "
                                   "give a negative argument for a logarithm.");
  }
  Real G0os = robust::ratio(_G0, _s);
  // sk, lk, mk, sum, and enough values may change in the loop
  Real sk = -_G0 * eta;
  Real lk = _G0;
  Real mk = 0.0;
  Real sum = sk;
  int enough = 0;
  Real temp;
  Real pf = _Cs * _Cs / (2.0 * _Cv0 * _s * _s);
  if (eta * eta < 1.e-8) {
    temp = _T0 * exp(_G0 * eta) +
           pf * (2.0 * std::log(f1) - _s * eta / (f1 * f1) * (3.0 * _s * eta - 2.0));
  } else {
    for (int i = 0; ((i < sumkmax) && (enough == 0)); i++) {
      mk = G0os / (i + 2) / (i + 2);
      sk = (sk * f1 - lk * eta) * (i + 1) * mk;
      lk = mk * (i + 1) * lk;
      sum = sum + sk;
      if (robust::ratio((sk * sk), (sum * sum)) < (cutoff * cutoff)) {
        enough = 1;
      }
    }
    if (enough == 0) {
#ifndef NDEBUG
      PORTABLE_WARN("Hugoniot Temperature not converged");
#endif // NDEBUG
    }
    // printf("sum=%e\n",sum);
    temp = _T0 - pf * ((G0os - 3.0) + exp(-G0os) * (G0os * G0os - 4.0 * G0os + 2.0) *
                                          (std::log(f1) + sum));
    temp = temp * exp(_G0 * eta);
    temp = temp - pf / f1 * ((4.0 - G0os) - 1.0 / f1);
  }
  return temp;
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real MGUsup::InternalEnergyFromDensityTemperature(
    const Real rho, const Real temp, Indexer_t &&lambda) const {
  Real value =
      HugInternalEnergyFromDensity(rho) + _Cv0 * (temp - HugTemperatureFromDensity(rho));
  return value;
}
template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real MGUsup::PressureFromDensityTemperature(
    const Real rho, const Real temp, Indexer_t &&lambda) const {
  Real value = HugPressureFromDensity(rho) +
               _G0 * _rho0 * _Cv0 * (temp - HugTemperatureFromDensity(rho));
  return value;
}
template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real MGUsup::EntropyFromDensityTemperature(
    const Real rho, const Real temp, Indexer_t &&lambda) const {
  Real eta = 1.0 - robust::ratio(_rho0, rho);
  Real value = _S0 - _G0 * _Cv0 * eta + _Cv0 * std::log(temp / _T0);
  return value;
}
template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real MGUsup::TExpansionCoeffFromDensityTemperature(
    const Real rho, const Real temp, Indexer_t &&lambda) const {
  Real value =
      robust::ratio(_Cv0 * _rho0 * _G0, TBulkModulusFromDensityTemperature(rho, temp));
  return value;
}
template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real MGUsup::TBulkModulusFromDensityTemperature(
    const Real rho, const Real temp, Indexer_t &&lambda) const {
  Real eta = 1.0 - robust::ratio(_rho0, rho);
  Real value = robust::ratio((1.0 + _s * eta - _G0 * _s * eta * eta), (1.0 - _s * eta));
  if (eta == 0.0) {
    value = _Cs * _Cs * _rho0;
  } else {
    value = value * robust::ratio(HugPressureFromDensity(rho), eta);
  }
  value = value - _G0 * _G0 * _Cv0 * _rho0 * HugTemperatureFromDensity(rho);
  value = robust::ratio(_rho0, rho) * value;
  return value;
}
template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real MGUsup::BulkModulusFromDensityTemperature(
    const Real rho, const Real temp, Indexer_t &&lambda) const {
  Real eta = 1.0 - robust::ratio(_rho0, rho);
  Real value =
      robust::ratio((1.0 + _s * eta - _G0 * _s * eta * eta), eta * (1.0 - _s * eta));
  if (eta == 0.0) {
    value = _Cs * _Cs * _rho0;
  } else {
    value = value * robust::ratio(HugPressureFromDensity(rho), eta);
  }
  value = value - _G0 * _G0 * _Cv0 * _rho0 * (temp - HugTemperatureFromDensity(rho));
  value = robust::ratio(_rho0, rho) * value;
  return value;
}
template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real MGUsup::TemperatureFromDensityInternalEnergy(
    const Real rho, const Real sie, Indexer_t &&lambda) const {
  Real value =
      (sie - HugInternalEnergyFromDensity(rho)) / _Cv0 + HugTemperatureFromDensity(rho);
  if (value < 0.0) {
#ifndef NDEBUG
    PORTABLE_WARN("Negative temperature");
#endif // NDEBUG
    value = 1.e-12;
  }
  return value;
}
template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real MGUsup::PressureFromDensityInternalEnergy(
    const Real rho, const Real sie, Indexer_t &&lambda) const {
  Real value = HugPressureFromDensity(rho) +
               _rho0 * _G0 * (sie - HugInternalEnergyFromDensity(rho));
  return value;
}
template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real
MGUsup::MinInternalEnergyFromDensity(const Real rho, Indexer_t &&lambda) const {
  MinInternalEnergyIsNotEnabled("MGUsup");
  return 0.0;
}
template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real MGUsup::EntropyFromDensityInternalEnergy(
    const Real rho, const Real sie, Indexer_t &&lambda) const {
  Real eta = 1.0 - robust::ratio(_rho0, rho);
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
template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real MGUsup::BulkModulusFromDensityInternalEnergy(
    const Real rho, const Real sie, Indexer_t &&lambda) const {
  Real eta = 1.0 - robust::ratio(_rho0, rho);
  Real value = robust::ratio((1.0 + _s * eta - _G0 * _s * eta * eta), (1.0 - _s * eta));
  if (eta == 0.0) {
    value = _Cs * _Cs * _rho0;
  } else {
    value = value * robust::ratio(HugPressureFromDensity(rho), eta);
  }
  value = value + _G0 * _G0 * _rho0 * (sie - HugInternalEnergyFromDensity(rho));
  value = robust::ratio(_rho0, rho) * value;
  return value;
}
// AEM: Give error since function is not well defined
template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION void MGUsup::DensityEnergyFromPressureTemperature(
    const Real press, const Real temp, Indexer_t &&lambda, Real &rho, Real &sie) const {
  EOS_ERROR("MGUsup::DensityEnergyFromPressureTemperature: "
            "Not implemented.\n");
}
// AEM: We should add entropy and Gruneissen parameters here so that it is complete
// If we add also alpha and BT, those should also be in here.
template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION void
MGUsup::FillEos(Real &rho, Real &temp, Real &sie, Real &press, Real &cv, Real &bmod,
                const unsigned long output, Indexer_t &&lambda) const {
  const unsigned long input = ~output; // everything that is not output is input
  if (thermalqs::density & output) {
    EOS_ERROR("MGUsup FillEos: Density is required input.\n");
  } else if (thermalqs::temperature & output &&
             thermalqs::specific_internal_energy & output) {
    EOS_ERROR("MGUsup FillEos: Density and Internal Energy or Density and Temperature "
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
template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION void
MGUsup::ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                               Real &bmod, Real &dpde, Real &dvdt,
                               Indexer_t &&lambda) const {
  // AEM: Added all variables I think should be output eventually
  Real tbmod;
  // Real entropy, alpha, Gamma;

  rho = _rho0;
  temp = _T0;
  sie = _E0;
  press = PressureFromDensityTemperature(rho, temp, lambda);
  // entropy = _S0;
  cv = _Cv0;
  tbmod = _Cs * _Cs * _rho0 - _G0 * _G0 * _Cv0 * _rho0 * _T0;
  // alpha = _A0;
  bmod = _Cs * _Cs * _rho0;
  // Gamma = _G0;
  // AEM: I suggest taking the two following away.
  dpde = _G0 * _rho0;
  dvdt = robust::ratio(-1.0, _T0 * _G0 * _rho0);
}
} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_MGUSUP_HPP_

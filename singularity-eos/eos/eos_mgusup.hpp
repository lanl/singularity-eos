//------------------------------------------------------------------------------
// © 2021-2025. Triad National Security, LLC. All rights reserved.  This
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
#include <singularity-eos/base/eos_error.hpp>
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
         const Real Cv0, const Real E0, const Real S0,
         const MeanAtomicProperties &AZbar = MeanAtomicProperties())
      : _rho0(rho0), _T0(T0), _Cs(Cs), _s(s), _G0(G0), _Cv0(Cv0), _E0(E0), _S0(S0),
        _AZbar(AZbar) {
    CheckParams();
  }
  PORTABLE_INLINE_FUNCTION void CheckParams() const;

  MGUsup GetOnDevice() { return *this; }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temp,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real
  PressureFromDensityTemperature(const Real rho, const Real temp,
                                 Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real
  EntropyFromDensityTemperature(const Real rho, const Real temp,
                                Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real EntropyFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityTemperature(
      const Real rho, const Real temp,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return _Cv0;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return _Cv0;
  }
  // added for testing AEM Dec 2023
  PORTABLE_INLINE_FUNCTION Real HugPressureFromDensity(const Real rho) const;
  PORTABLE_INLINE_FUNCTION Real HugInternalEnergyFromDensity(const Real rho) const;
  PORTABLE_INLINE_FUNCTION Real HugTemperatureFromDensity(const Real rho) const;
  // Thermal Bulk Modulus added AEM Dec 2022
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real TBulkModulusFromDensityTemperature(
      const Real rho, const Real temp,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityTemperature(
      const Real rho, const Real temp,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  // Thermal expansion coefficient added AEM 2022
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real TExpansionCoeffFromDensityTemperature(
      const Real rho, const Real temp,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temp,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return robust::ratio(_G0 * _rho0, rho);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return robust::ratio(_G0 * _rho0, rho);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  FillEos(Real &rho, Real &temp, Real &energy, Real &press, Real &cv, Real &bmod,
          const unsigned long output,
          Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;

  // Since reference isotherm scales as eta = 1 -_rho0/rho, EOS
  // diverges for rho << rho0. eta^2 is used, so...
  PORTABLE_FORCEINLINE_FUNCTION
  Real MinimumDensity() const { return _RHOMINFAC * _rho0; }

  PORTABLE_FORCEINLINE_FUNCTION
  Real MaximumDensity() const {
    if (_s > 1) {
      return 0.99 * robust::ratio(_s * _rho0, _s - 1);
    } else { // for s <= 1, no maximum, but we need to pick something.
      return 1e3 * _rho0;
    } // note that s < 0 implies unphysical shock derivative.
  }
  // Hugoniot pressure ill defined at reference density. On one side,
  // negative. On the other positive.
  PORTABLE_FORCEINLINE_FUNCTION
  Real MinimumPressure() const { return -1e100; }
  PORTABLE_FORCEINLINE_FUNCTION
  Real MaximumPressureAtTemperature([[maybe_unused]] const Real T) const { return 1e100; }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                         Real &bmod, Real &dpde, Real &dvdt,
                         Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  // Generic functions provided by the base class. These contain e.g. the vector
  // overloads that use the scalar versions declared here
  SG_ADD_DEFAULT_MEAN_ATOMIC_FUNCTIONS(_AZbar)
  SG_ADD_BASE_CLASS_USINGS(MGUsup)

  static constexpr unsigned long PreferredInput() { return _preferred_input; }
  PORTABLE_INLINE_FUNCTION void PrintParams() const {
    static constexpr char st[]{"MGUsup Params: "};
    printf("%s rho0:%e T0:%e Cs:%e s:%e\n  G0:%e Cv0:%e E0:%e S0:%e\n", st, _rho0, _T0,
           _Cs, _s, _G0, _Cv0, _E0, _S0);
    _AZbar.PrintParams();
    printf("\n\n");
  }
  inline void Finalize() {}
  static std::string EosType() { return std::string("MGUsup"); }
  static std::string EosPyType() { return EosType(); }

 private:
  static constexpr const unsigned long _preferred_input =
      thermalqs::density | thermalqs::specific_internal_energy;
  Real _RHOMINFAC = std::sqrt(robust::EPS());
  Real _rho0, _T0, _Cs, _s, _G0, _Cv0, _E0, _S0;
  MeanAtomicProperties _AZbar;
};

PORTABLE_INLINE_FUNCTION void MGUsup::CheckParams() const {
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
  _AZbar.CheckParams();
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
    printf("f1, eta, rho, rho0, s = %.14e %.14e %.14e %.14e %.14e\n", f1, eta, rho, _rho0,
           _s);
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
  // if (output & thermalqs::specific_internal_energy) sie = sie;
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
  // Real tbmod;
  // Real entropy, alpha, Gamma;

  rho = _rho0;
  temp = _T0;
  sie = _E0;
  press = PressureFromDensityTemperature(rho, temp, lambda);
  // entropy = _S0;
  cv = _Cv0;
  // tbmod = _Cs * _Cs * _rho0 - _G0 * _G0 * _Cv0 * _rho0 * _T0;
  // alpha = _A0;
  bmod = _Cs * _Cs * _rho0;
  // Gamma = _G0;
  // AEM: I suggest taking the two following away.
  dpde = _G0 * _rho0;
  dvdt = robust::ratio(-1.0, _T0 * _G0 * _rho0);
}

#ifdef SINGULARITY_INSTANTIATE_CLASSES
SG_ADD_TEMPLATE_EXTERNS(MGUsup, Real *)
#endif // SINGULARITY_INSTANTIATE_CLASSES

} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_MGUSUP_HPP_

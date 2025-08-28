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

#ifndef _SINGULARITY_EOS_EOS_EOS_SIMPLE_MACAW_HPP_
#define _SINGULARITY_EOS_EOS_EOS_SIMPLE_MACAW_HPP_

// stdlib
#include <cmath>
#include <cstdio>
#include <string>

// Ports-of-call
#include <ports-of-call/portability.hpp>

// Base stuff
#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/eos_error.hpp>
#include <singularity-eos/base/robust_utils.hpp>
#include <singularity-eos/eos/eos_base.hpp>

namespace singularity {

using namespace eos_base;

class SimpleMACAW : public EosBase<SimpleMACAW> {
 public:
  SimpleMACAW() = default;
  PORTABLE_INLINE_FUNCTION
  SimpleMACAW(Real A, Real B, Real Cvinf, Real v0, Real T0, Real Gc,
           const MeanAtomicProperties &AZbar = MeanAtomicProperties())
      : _A(A), _B(B), _Cvinf(Cvinf), _v0(v0), _T0(T0), _Gc(Gc), _AZbar(AZbar) {
    CheckParams();
  }
  SimpleMACAW GetOnDevice() { return *this; }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real x = std::pow( rho * _v0, _Gc);
    const Real sie_term = (sie - SieColdCurve(1.0 / rho)) / _Cvinf;
    // It's a quadratic equation to solve for temperature from sie and rho
    const Real b = -sie_term;
    const Real c = -x * sie_term;
    return 0.5 * (-b + std::sqrt(b * b - 4.0 * c));
  }
  PORTABLE_INLINE_FUNCTION void CheckParams() const {
    PORTABLE_ALWAYS_REQUIRE(_A > 0, "Parameter, 'A', must be positive");
    PORTABLE_ALWAYS_REQUIRE(_B > 0, "Parameter, 'B', must be positive");
    PORTABLE_ALWAYS_REQUIRE(_v0 > 0, "Reference specific volume, 'v0', must be positive");
    PORTABLE_ALWAYS_REQUIRE(_T0 > 0, "Reference temperature, 'T0', must be positive");
    PORTABLE_ALWAYS_REQUIRE(_Cvinf > 0, "Specific heat capacity (at constant volume), `Cvinf`, must be positive");
    PORTABLE_ALWAYS_REQUIRE(_Gc > 0 && _Gc < 1, "Gruneisen coefficient, 'Gc', must be in the interval (0,1)");
    _AZbar.CheckParams();
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real v = 1.0 / rho;
    return rho * (SieColdCurve(v) + SieThermalPortion(v, temperature));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real v = 1.0 / rho;
    return PressureColdCurve(v) + PressureThermalPortion(v, temperature);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const v = 1.0 / rho;
    return PressureColdCurve(v) + _Gc * rho * (sie - SieColdCurve(v));
  }

  // Entropy member functions
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real
  EntropyFromDensityTemperature(const Real rho, const Real temperature,
                                Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real tau = robust::ratio(temperature, TemperatureScale(rho));
    return _Cvinf * (robust::ratio(tau, 1.0 + tau) + std::log(1.0 + tau));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real EntropyFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real T = TemperatureFromDensityInternalEnergy(rho, sie);
    return EntropyFromDensityTemperature(rho, T);
  }

  // Cold curve, thermal portion, and other helper functions
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SieColdCurve(
      const Real v, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real x = robust::ratio(v, _v0);
    return _A * _v0 * (std::pow(x, -_B) + _B * x - (_B + 1.0));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SieThermalPortion(
      const Real v, const Real T,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real x = robust::ratio(v, _v0);
    return _Cvinf * T * (1.0 - _T0 / (_T0 + T * std::pow(x, _Gc)));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureColdCurve(
      const Real v, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real x = robust::ratio(v, _v0);
    return _A * _B * (std::pow(x, -(_B + 1.0)) - 1.0);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureThermalPortion(
      const Real v, const Real T,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return rho * _Cvinf * _Gc * pow<2>(T) / (_T0 * std::pow(rho * _v0, _Gc) + T);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real TemperatureScale(
      const Real rho, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return _T0 * std::pow(rho * _v0, _Gc);
  }

  // Thermodynamic derivatives
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real tau = robust::ratio(temperature, TemperatureScale(rho));
    const Real numerator = pow<2>(tau) + 2.0 * tau;
    const Real denominator = pow<2>(tau + 1.0);
    return _Cvinf * robust::ratio(numerator / denominator);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real T = TemperatureFromDensityInternalEnergy(rho, sie);
    return SpecificHeatFromDensityTemperature(rho, T);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return std::max(robust::SMALL(), _gm1 * (_gm1 + 1.0) * rho * _Cv * temperature);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real term1 = _A * _B * (_B + 1.0) * std::pow(rho * _v0, _B + 1.0) 
                       + _Gc * (_Gc + 1.0) * rho * (sie - SieColdCurve(1.0 / rho));
    return std::max(robust::SMALL(), _gm1 * (_gm1 + 1.0) * (rho * (sie - _qq) - _Pinf));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return _Gc;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return _Gc;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  FillEos(Real &rho, Real &temp, Real &energy, Real &press, Real &cv, Real &bmod,
          const unsigned long output,
          Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  //template <typename Indexer_t = Real *>
  //PORTABLE_INLINE_FUNCTION void
  //ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
  //                       Real &bmod, Real &dpde, Real &dvdt,
  //                       Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
  //  // use STP: 1 atmosphere, room temperature
  //  rho = _rho0;
  //  temp = _T0;
  //  sie = _sie0;
  //  press = _P0;
  //  cv = _Cv;
  //  bmod = _bmod0;
  //  dpde = _dpde0;
  //}
  // Generic functions provided by the base class. These contain e.g. the vector
  // overloads that use the scalar versions declared here
  SG_ADD_BASE_CLASS_USINGS(SimpleMACAW)
  SG_ADD_DEFAULT_MEAN_ATOMIC_FUNCTIONS(_AZbar)

  static constexpr unsigned long PreferredInput() { return _preferred_input; }
  PORTABLE_INLINE_FUNCTION void PrintParams() const {
    printf("Simple MACAW Parameters:\nA = %g\nB    = %g\nCvinf  = %g\nGc     = "
           "%g\n",
           _A, _B, _Cvinf, _Gc);
    _AZbar.PrintParams();
  }
  template <typename Indexer_t>
  PORTABLE_INLINE_FUNCTION void
  DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                       Indexer_t &&lambda, Real &rho, Real &sie) const {
    sie = std::max(
        _qq,
        robust::ratio(press + (_gm1 + 1.0) * _Pinf, press + _Pinf) * _Cv * temp + _qq);
    rho = std::max(robust::SMALL(), robust::ratio(press + _Pinf, _gm1 * _Cv * temp));
  }
  inline void Finalize() {}
  static std::string EosType() { return std::string("SimpleMACAW"); }
  static std::string EosPyType() { return EosType(); }

 private:
  Real _A, _B, _Cvinf, _Gc;
  // reference values
  Real _v0, _T0;
  MeanAtomicProperties _AZbar;
  static constexpr const unsigned long _preferred_input =
      thermalqs::density | thermalqs::specific_internal_energy;
};

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION void
StiffGas::FillEos(Real &rho, Real &temp, Real &sie, Real &press, Real &cv, Real &bmod,
                  const unsigned long output, Indexer_t &&lambda) const {
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
    sie = robust::ratio(press + (_gm1 + 1.0) * _Pinf, _gm1 * rho) + _qq;
  }
  if (output & thermalqs::pressure) press = PressureFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::temperature)
    temp = TemperatureFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::bulk_modulus)
    bmod = BulkModulusFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::specific_heat)
    cv = SpecificHeatFromDensityInternalEnergy(rho, sie);
}

} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_EOS_SIMPLE_MACAW_HPP_

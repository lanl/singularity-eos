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

/* Most information regarding this equation of state can be found here:
 * https://www.osti.gov/biblio/2479474 
 * 
 * Note that all equations are given as functions of the density ($\rho$) rather than functions
 * of the specific volume ($v$). The reason for this is it simplifies the computational complexity.
 */
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
    /* Equation (18) */
    const Real Delta_e = sie - SieColdCurve(rho);
    const Real discriminant = Delta_e * (Delta_e + 4.0 * _Cvinf * _T0 * std::pow(rho * _v0, _Gc));
    return robust::ratio((Delta_e + std::sqrt(discriminant)), 2.0 * _Cvinf);
  }
  PORTABLE_INLINE_FUNCTION void CheckParams() const {
    PORTABLE_ALWAYS_REQUIRE(_A >= 0, "Parameter, 'A', must be non-negative");
    PORTABLE_ALWAYS_REQUIRE(_B >= 0, "Parameter, 'B', must be non-negative");
    PORTABLE_ALWAYS_REQUIRE(_Cvinf > 0, "Specific heat capacity (at constant volume), `Cvinf`, must be positive");
    PORTABLE_ALWAYS_REQUIRE(_v0 > 0, "Reference specific volume, 'v0', must be positive");
    PORTABLE_ALWAYS_REQUIRE(_T0 >= 0, "Reference temperature, 'T0', must be non-negative");
    PORTABLE_ALWAYS_REQUIRE(_Gc > 0, "Gruneisen parameter, 'Gc', must be positive");
    if (_Gc > 1) {PORTABLE_WARN("Warning: Gruneisen coefficient greater than 1. Thermodynamic stability may be violated.");}
    _AZbar.CheckParams();
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
      /* Equation (16) */
    return rho * (SieColdCurve(rho) + SieThermalPortion(rho, temperature));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
      /* Equation (15) */
    return PressureColdCurve(rho) + PressureThermalPortion(rho, temperature);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
      /* Equation (19) */
    return PressureColdCurve(rho) + _Gc * rho * (sie - SieColdCurve(rho));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real InternalEnergyFromDensityPressure(
      const Real rho, const Real P,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return SieColdCurve(rho) + robust::ratio(P - PressureColdCurve(rho), _Gc * rho);
  }

  // Entropy member functions
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real
  EntropyFromDensityTemperature(const Real rho, const Real temperature,
                                Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    /* Equation (8) */
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

  /* Note: The cold curves for the simple MACAW are presented as functions of `v` rather than `rho`.
   * We keep them as functions of `rho` for computational efficiency. */
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SieColdCurve(
      const Real rho, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real ratio = rho * _v0;
    return _A * _v0 * (std::pow(ratio, _B) + robust::ratio(_B, ratio) - (_B + 1.));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SieThermalPortion(
      const Real rho, const Real T,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real ratio = rho * _v0;
    return _Cvinf * T * (1.0 - robust::ratio(_T0, _T0 + T * std::pow(ratio, -_Gc)));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureColdCurve(
      const Real rho, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real ratio = rho * _v0;
    return _A * _B * (std::pow(ratio, _B + 1.0) - 1.0);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureThermalPortion(
      const Real rho, const Real T,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
      const Real numerator = rho * _Cvinf * _Gc * math_utils::pow<2>(T);
      const Real denominator = _T0 * std::pow(rho * _v0, _Gc) + T;
    return robust::ratio(numerator, denominator);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real TemperatureScale(
      const Real rho, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    /* Equation (13) */
    return _T0 * std::pow(rho * _v0, _Gc);
  }

  // Specific heat capacity at constant volume
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    /* Equation (6) */
    const Real tau = robust::ratio(temperature, TemperatureScale(rho));
    const Real numerator = math_utils::pow<2>(tau) + 2.0 * tau;
    const Real denominator = math_utils::pow<2>(tau + 1.0);
    return _Cvinf * robust::ratio(numerator, denominator);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real T = TemperatureFromDensityInternalEnergy(rho, sie);
    return SpecificHeatFromDensityTemperature(rho, T);
  }
  // Isentropic bulk modulus 
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real sie = InternalEnergyFromDensityTemperature(rho, temperature);
    return BulkModulusFromDensityInternalEnergy(rho, sie);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    /* Equation (21) (multiplied by $\rho$ to get bulk mod) */
    const Real term1 = _A * _B * (_B + 1.0) * std::pow(rho * _v0, _B + 1.0);
    const Real term2 = _Gc * (_Gc + 1.0) * rho * (sie - SieColdCurve(rho));
    return term1 + term2;
  }
  // Isothermal bulk modulus
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real IsothermalBulkModulusFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    /* Equation (22) (multiplied by $\rho$ to get bulk mod)
     * As mentioned in the paper, from the constraints on the parameters,
     * the isothermal bulk modulus is guaranteed to be positive. */
    const Real term1 = _A * _B * (_B + 1.0) * std::pow(rho * _v0, _B + 1.0);
    const Real numerator = rho * _T0 * (1.0 - _Gc) * std::pow(rho * _v0, 2.0 * _Gc) + temperature;
    const Real denominator = _T0 * std::pow(rho * _v0, _Gc) + temperature;
    return term1 + _Cvinf * _Gc * temperature * temperature * robust::ratio(numerator, denominator * denominator);
  }
  // Specific heat capacity at constant pressure
   template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real ConstantPressureSpecificHeatFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real BT = IsothermalBulkModulusFromDensityTemperature(rho, temperature);
    const Real Bs = BulkModulusFromDensityTemperature(rho, temperature);
    const Real cv = SpecificHeatFromDensityTemperature(rho, temperature);
    return robust::ratio(Bs * cv, BT); /* General thermodynamic identity */
  }

  // Gruneisen parameter
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
  FillEos(Real &rho, Real &temp, Real &sie, Real &press, Real &cv, Real &bmod,
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
    sie = InternalEnergyFromDensityPressure(rho, press);
  }
  if (output & thermalqs::pressure) press = PressureFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::temperature)
    temp = TemperatureFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::bulk_modulus)
    bmod = BulkModulusFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::specific_heat)
    cv = SpecificHeatFromDensityInternalEnergy(rho, sie);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                         Real &bmod, Real &dpde, Real &dvdt,
                         Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    // v_0 and T_0 are the only user selected reference state parameters
    rho = 1.0 / _v0;
    temp = _T0;
    press = PressureFromDensityTemperature(rho, temp);
    sie = InternalEnergyFromDensityTemperature(rho, temp);

    cv = SpecificHeatFromDensityInternalEnergy(rho, sie);
    bmod = BulkModulusFromDensityInternalEnergy(rho, sie);
    dpde = rho * GruneisenParamFromDensityTemperature(rho, temp);
  }
  // Generic functions provided by the base class. These contain e.g. the vector
  // overloads that use the scalar versions declared here
  SG_ADD_BASE_CLASS_USINGS(SimpleMACAW)
  SG_ADD_DEFAULT_MEAN_ATOMIC_FUNCTIONS(_AZbar)

  static constexpr unsigned long PreferredInput() { return _preferred_input; }
  PORTABLE_INLINE_FUNCTION void PrintParams() const {
    printf("Simple MACAW Parameters:\nA = %g\nB    = %g\nCvinf  = %g\nv0  = %g\nT0  = %g\nGc     = "
           "%g\n",
           _A, _B, _Cvinf, _v0, _T0, _Gc);
    _AZbar.PrintParams();
  }
  template <typename Indexer_t>
  PORTABLE_INLINE_FUNCTION void
  DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                       Indexer_t &&lambda, Real &rho, Real &sie) const {
    /* Since the isothermal bulk modulus is always positive (assuming parameter constraints),
     * this guarantees the function to solve is monotonic hence we always have a unique solution. */

    // Setup lambda function for rho
    auto f = [=](const Real x /* density */){
      const Real term1 = _A * _B * (std::pow(_v0 * x, _B + 1.0) - 1.0) - press;
      const Real term2 = _T0 * std::pow(_v0 * x, _Gc) + temp;
      const Real term3 = _Cvinf * _Gc * temp * temp * x;
      return term1 * term2 + term3;
    };

    const RootFinding1D::RootCounts root_info;

    // Finding the density on pressure cold curve is a guaranteed upper bound
    const Real rho_high = std::pow(press / (_A * _B) + 1.0, 1.0 / (_B + 1.0)) / _v0;
    const Real rho_low = 0.0; // zero density is valid for `f` defined above.

    regula_falsi(f, 0.0 /*target*/, 1.0 / _v0 /*guess*/, 
                 rho_low /*left bracket*/, rho_high /*right bracket*/,
                 1.0e-9 /*x? tol*/, 1.0e-9 /*y? tol*/,
                 rho, &root_info, true);

    sie = InternalEnergyFromDensityTemperature(rho, temp);
  }
  inline void Finalize() {}
  static std::string EosType() { return std::string("SimpleMACAW"); }
  static std::string EosPyType() { return EosType(); }

 private:
  Real _A, _B, _Cvinf, _v0, _T0, _Gc;
  MeanAtomicProperties _AZbar;
  static constexpr const unsigned long _preferred_input =
      thermalqs::density | thermalqs::specific_internal_energy;
};

} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_EOS_SIMPLE_MACAW_HPP_

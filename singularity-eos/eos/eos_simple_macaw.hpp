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
 * Note that all equations are given as functions of the density ($\rho$) rather than
 * functions of the specific volume ($v$). The reason for this is it simplifies the
 * computational complexity.
 */
class SimpleMACAW : public EosBase<SimpleMACAW> {
 public:
  SimpleMACAW() = default;
  PORTABLE_INLINE_FUNCTION
  SimpleMACAW(Real A, Real B, Real Cvinf, Real v0, Real T0, Real Gc,
              const MeanAtomicProperties &AZbar = MeanAtomicProperties())
      : A_(A), B_(B), Cvinf_(Cvinf), v0_(v0), T0_(T0), Gc_(Gc), _AZbar(AZbar) {
    CheckParams();
  }
  SimpleMACAW GetOnDevice() { return *this; }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    /* Equation (18) */
    const Real Delta_e = std::max(sie - SieColdCurve(rho), 0.);

    // Early return ensure we're on the cold curve
    if (_IsNearOrBelowZero(Delta_e)) {
      return 0.;
    }

    const Real discriminant =
        Delta_e * (Delta_e + 4.0 * Cvinf_ * T0_ * std::pow(rho * v0_, Gc_));

    return robust::ratio((Delta_e + std::sqrt(discriminant)), 2.0 * Cvinf_);
  }
  PORTABLE_INLINE_FUNCTION void CheckParams() const {
    PORTABLE_ALWAYS_REQUIRE(A_ >= 0, "Parameter, 'A', must be non-negative");
    PORTABLE_ALWAYS_REQUIRE(B_ >= 0, "Parameter, 'B', must be non-negative");
    PORTABLE_ALWAYS_REQUIRE(
        Cvinf_ > 0,
        "Specific heat capacity (at constant volume), `Cvinf`, must be positive");
    PORTABLE_ALWAYS_REQUIRE(v0_ > 0, "Reference specific volume, 'v0', must be positive");
    PORTABLE_ALWAYS_REQUIRE(T0_ >= 0,
                            "Reference temperature, 'T0', must be non-negative");
    PORTABLE_ALWAYS_REQUIRE(Gc_ > 0, "Gruneisen parameter, 'Gc', must be positive");
    if (Gc_ > 1) {
      PORTABLE_WARN("Warning: Gruneisen coefficient greater than 1. Thermodynamic "
                    "stability may be violated.");
    }
    _AZbar.CheckParams();
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    /* Equation (16) */
    return SieColdCurve(rho) + _SieThermalPortion(rho, temperature);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    /* Equation (15) */
    return PressureColdCurve(rho) + _PressureThermalPortion(rho, temperature);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    /* Equation (19) */
    const Real Delta_e = std::max(sie - SieColdCurve(rho), 0.);

    // Early return ensure we're on the cold curve
    if (_IsNearOrBelowZero(Delta_e)) {
      return PressureColdCurve(rho);
    }

    return PressureColdCurve(rho) + Gc_ * rho * Delta_e;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real InternalEnergyFromDensityPressure(
      const Real rho, const Real P,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real delta_p = std::max(P - PressureColdCurve(rho), 0.);

    // Early return ensure we're on the cold curve
    if (_IsNearOrBelowZero(delta_p)) {
      return SieColdCurve(rho);
    }

    return SieColdCurve(rho) + robust::ratio(delta_p, Gc_ * rho);
  }

  // Entropy member functions
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real
  EntropyFromDensityTemperature(const Real rho, const Real temperature,
                                Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    /* Equation (8) */
    const Real tau = robust::ratio(temperature, _TemperatureScale(rho));
    return Cvinf_ * (robust::ratio(tau, 1.0 + tau) + std::log(1.0 + tau));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real EntropyFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real T = TemperatureFromDensityInternalEnergy(rho, sie);
    return EntropyFromDensityTemperature(rho, T);
  }

  // Cold curve, thermal portion, and other helper functions

  /* Note: The cold curves for the simple MACAW are presented as functions of `v` rather
   * than `rho`. We keep them as functions of `rho` for computational efficiency. */
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real
  SieColdCurve(const Real rho, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real ratio = rho * v0_;
    return A_ * v0_ * (std::pow(ratio, B_) + robust::ratio(B_, ratio) - (B_ + 1.));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureColdCurve(
      const Real rho, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real ratio = rho * v0_;
    return A_ * B_ * (std::pow(ratio, B_ + 1.0) - 1.0);
  }

  // Specific heat capacity at constant volume
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    /* Equation (6) */
    const Real tau = robust::ratio(temperature, _TemperatureScale(rho));
    const Real numerator = math_utils::pow<2>(tau) + 2.0 * tau;
    const Real denominator = math_utils::pow<2>(tau + 1.0);
    return Cvinf_ * robust::ratio(numerator, denominator);
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
    const Real Delta_e = std::max(sie - SieColdCurve(rho), 0.);
    /* Equation (21) (multiplied by $\rho$ to get bulk mod) */
    const Real term1 = A_ * B_ * (B_ + 1.0) * std::pow(rho * v0_, B_ + 1.0);

    // Early return ensure we're on the cold curve
    if (_IsNearOrBelowZero(Delta_e)) {
      return term1;
    }

    const Real term2 = Gc_ * (Gc_ + 1.0) * rho * Delta_e;
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
    const Real term1 = A_ * B_ * (B_ + 1.0) * std::pow(rho * v0_, B_ + 1.0);
    // There is a mistake in the paper and the numerator term should have a -Gc_ instead
    // of -2*Gc_
    const Real numerator =
        rho * (T0_ * (1.0 - Gc_) * std::pow(rho * v0_, Gc_) + temperature);
    const Real denominator = T0_ * std::pow(rho * v0_, Gc_) + temperature;
    return term1 + Cvinf_ * Gc_ * temperature * temperature *
                       robust::ratio(numerator, denominator * denominator);
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
  // Coefficient of thermal expansivity
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real CoefficientThermalExpansionFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real BT = IsothermalBulkModulusFromDensityTemperature(rho, temperature);
    const Real cv = SpecificHeatFromDensityTemperature(rho, temperature);
    return robust::ratio(Gc_ * rho * cv, BT); /* General thermodynamic identity */
  }

  // Gruneisen parameter
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return Gc_;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return Gc_;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void FillEos(Real &rho, Real &temp, Real &sie, Real &press,
                                        Real &cv, Real &bmod, const unsigned long output,
                                        Indexer_t &&lambda) const {
    if (output & thermalqs::density && output & thermalqs::specific_internal_energy) {
      PORTABLE_ALWAYS_REQUIRE(
          !(output & thermalqs::pressure || output & thermalqs::temperature),
          "Cannot request more than two thermodynamic variables (rho, T, e, P)");
      DensityEnergyFromPressureTemperature(press, temp, lambda, rho, sie);
    }
    if (output & thermalqs::pressure && output & thermalqs::specific_internal_energy) {
      PORTABLE_ALWAYS_REQUIRE(
          !(output & thermalqs::density || output & thermalqs::temperature),
          "Cannot request more than two thermodynamic variables (rho, T, e, P)");
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
    rho = 1.0 / v0_;
    temp = T0_;
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
    printf("Simple MACAW Parameters:\nA = %g\nB    = %g\nCvinf  = %g\nv0  = %g\nT0  = "
           "%g\nGc     = "
           "%g\n",
           A_, B_, Cvinf_, v0_, T0_, Gc_);
    _AZbar.PrintParams();
  }
  template <typename Indexer_t>
  PORTABLE_INLINE_FUNCTION void
  DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                       Indexer_t &&lambda, Real &rho, Real &sie) const {
    /* Since the isothermal bulk modulus is always positive (assuming parameter
     * constraints), this guarantees the function to solve is monotonic hence we always
     * have a unique solution. */

    // Setup lambda function for rho
    auto f = [=](const Real x /* density */) {
      const Real term1 = A_ * B_ * (std::pow(v0_ * x, B_ + 1.0) - 1.0) - press;
      const Real term2 = T0_ * std::pow(v0_ * x, Gc_) + temp;
      const Real term3 = Cvinf_ * Gc_ * temp * temp * x;
      return term1 * term2 + term3;
    };

    const RootFinding1D::RootCounts root_info;

    // Finding the density on pressure cold curve is a guaranteed upper bound
    const Real rho_high = std::pow(press / (A_ * B_) + 1.0, 1.0 / (B_ + 1.0)) / v0_;
    const Real rho_low = 0.0; // zero density is valid for `f` defined above.

    regula_falsi(f, 0.0 /*target*/, 1.0 / v0_ /*guess*/, rho_low /*left bracket*/,
                 rho_high /*right bracket*/, 1.0e-9 /*x? tol*/, 1.0e-9 /*y? tol*/, rho,
                 &root_info, true);

    sie = InternalEnergyFromDensityTemperature(rho, temp);
  }
  inline void Finalize() {}
  static std::string EosType() { return std::string("SimpleMACAW"); }
  static std::string EosPyType() { return EosType(); }

 private:
  Real A_, B_, Cvinf_, v0_, T0_, Gc_;
  MeanAtomicProperties _AZbar;
  static constexpr const unsigned long _preferred_input =
      thermalqs::density | thermalqs::specific_internal_energy;

  PORTABLE_FORCEINLINE_FUNCTION bool _IsNearOrBelowZero(Real value) {
    return value < std::numeric_limits<Real>min() * 5;
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real _TemperatureScale(
      const Real rho, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    /* Equation (13) */
    return T0_ * std::pow(rho * v0_, Gc_);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real
  _SieThermalPortion(const Real rho, const Real T,
                     Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real ratio = rho * v0_;
    return Cvinf_ * T * (1.0 - robust::ratio(T0_, T0_ + T * std::pow(ratio, -Gc_)));
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real
  _PressureThermalPortion(const Real rho, const Real T,
                          Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real numerator = rho * Cvinf_ * Gc_ * math_utils::pow<2>(T);
    const Real denominator = T0_ * std::pow(rho * v0_, Gc_) + T;
    return robust::ratio(numerator, denominator);
  }
};

} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_EOS_SIMPLE_MACAW_HPP_

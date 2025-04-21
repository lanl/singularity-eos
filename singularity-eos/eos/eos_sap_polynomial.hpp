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

#ifndef _SINGULARITY_EOS_EOS_EOS_SAP_POLYNOMIAL_HPP_
#define _SINGULARITY_EOS_EOS_EOS_SAP_POLYNOMIAL_HPP_

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

class SAP_Polynomial : public EosBase<SAP_Polynomial> {
 public:
  SAP_Polynomial() = default;
  PORTABLE_INLINE_FUNCTION
  SAP_Polynomial(const Real rho0, const Real a0, const Real a1, const Real a2c,
                 const Real a2e, const Real a3, const Real b0, const Real b1,
                 const Real b2c, const Real b2e, const Real b3,
                 const MeanAtomicProperties &AZbar = MeanAtomicProperties())
      : _a0(a0), _a1(a1), _a2c(a2c), _a2e(a2e), _a3(a3), _b0(b0), _b1(b1), _b2c(b2c),
        _b2e(b2e), _b3(b3), _rho0(rho0), _AZbar(AZbar) {
    CheckParams();
  }

  SAP_Polynomial GetOnDevice() { return *this; }
  PORTABLE_INLINE_FUNCTION void CheckParams() const {
    PORTABLE_ALWAYS_REQUIRE(_rho0 >= 0, "Reference density must be non-negative");
    _AZbar.CheckParams();
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    PORTABLE_WARN("This function is a stub for an incomplete EoS.");
    return 0.0;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    PORTABLE_WARN("This function is a stub for an incomplete EoS.");
    return 0.0;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    PORTABLE_WARN("This function is a stub for an incomplete EoS.");
    return 0.0;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real mu = MuFromDensity(rho);
    if (mu >= 0) // Compression
      return _a0 + _a1 * mu + _a2c * mu * mu + _a3 * mu * mu * mu +
             sie * (_b0 + _b1 * mu + _b2c * mu * mu + _b3 * mu * mu * mu);
    else
      return _a0 + _a1 * mu + _a2e * mu * mu + _a3 * mu * mu * mu +
             sie * (_b0 + _b1 * mu + _b2e * mu * mu + _b3 * mu * mu * mu);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real
  EntropyFromDensityTemperature(const Real rho, const Real temperature,
                                Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    PORTABLE_WARN("This function is a stub for an incomplete EoS.");
    return 0.0;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real EntropyFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    PORTABLE_WARN("This function is a stub for an incomplete EoS.");
    return 0.0;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    PORTABLE_WARN("This function is a stub for an incomplete EoS.");
    return 0.0;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    PORTABLE_WARN("This function is a stub for an incomplete EoS.");
    return 0.0;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    PORTABLE_WARN("This function is a stub for an incomplete EoS.");
    return 0.0;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    PORTABLE_WARN("This function is a stub for an incomplete EoS.");
    return 0.0;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real mu = MuFromDensity(rho);
    if (mu >= 0) // Compression
      return _b0 + _b1 * mu + _b2c * mu * mu + _b3 * mu * mu * mu;
    else
      return _b0 + _b1 * mu + _b2e * mu * mu + _b3 * mu * mu * mu;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real mu = MuFromDensity(rho);
    if (mu >= 0) // Compression
      return (1 + mu) * (_a1 + 2 * _a2c * mu + 3 * _a3 * mu * mu +
                         sie * (_b1 + 2 * _b2c * mu + 3 * _b3 * mu * mu));
    else
      return (1 + mu) * (_a1 + 2 * _a2e * mu + 3 * _a3 * mu * mu +
                         sie * (_b1 + 2 * _b2e * mu + 3 * _b3 * mu * mu));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real
  MuFromDensity(const Real rho, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return rho / _rho0 - 1;
  }

  // Essentially unbounded... I think.
  PORTABLE_FORCEINLINE_FUNCTION
  Real MinimumPressure() const { return -1e100; }
  PORTABLE_FORCEINLINE_FUNCTION
  Real MaximumPressureAtTemperature([[maybe_unused]] const Real T) const { return 1e100; }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  FillEos(Real &rho, Real &temp, Real &energy, Real &press, Real &cv, Real &bmod,
          const unsigned long output,
          Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                         Real &bmod, Real &dpde, Real &dvdt,
                         Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    // use STP: 1 atmosphere, room temperature
    rho = _rho0;
    temp = 0.0;
    sie = 0.0;
    press = 0.0;
    cv = 0.0;
    bmod = 0.0;
    dpde = 0.0;
    dvdt = 0.0;
  }
  // Generic functions provided by the base class. These contain e.g. the vector
  // overloads that use the scalar versions declared here
  SG_ADD_DEFAULT_MEAN_ATOMIC_FUNCTIONS(_AZbar)
  SG_ADD_BASE_CLASS_USINGS(SAP_Polynomial)

  static constexpr unsigned long PreferredInput() { return _preferred_input; }
  PORTABLE_INLINE_FUNCTION void PrintParams() const {
    printf("SAP_Polynomial EOS Parameters:\n");
    printf("    rho0  = %g\n", _rho0);
    printf("      a0  = %g\n", _a0);
    printf("      a1  = %g\n", _a1);
    printf("      a2c = %g\n", _a2c);
    printf("      a2e = %g\n", _a2e);
    printf("      a3  = %g\n", _a3);
    printf("      b0  = %g\n", _b0);
    printf("      b1  = %g\n", _b1);
    printf("      b2c = %g\n", _b2c);
    printf("      b2e = %g\n", _b2e);
    printf("      b3  = %g\n", _b3);
    _AZbar.PrintParams();
  }
  inline void Finalize() {}
  static std::string EosType() { return std::string("SAP_Polynomial"); }
  static std::string EosPyType() { return EosType(); }

 private:
  Real _a0, _a1, _a2c, _a2e, _a3, _b0, _b1, _b2c, _b2e, _b3;
  // reference values
  Real _rho0;
  MeanAtomicProperties _AZbar;
  static constexpr const unsigned long _preferred_input =
      thermalqs::density | thermalqs::specific_internal_energy;
};

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION void SAP_Polynomial::FillEos(Real &rho, Real &temp, Real &sie,
                                                      Real &press, Real &cv, Real &bmod,
                                                      const unsigned long output,
                                                      Indexer_t &&lambda) const {
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
    sie = 0.0;
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

#endif // _SINGULARITY_EOS_EOS_EOS_SAP_POLYNOMIAL_HPP_

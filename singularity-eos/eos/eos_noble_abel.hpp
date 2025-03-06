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

#ifndef _SINGULARITY_EOS_EOS_EOS_NOBLE_ABEL_HPP_
#define _SINGULARITY_EOS_EOS_EOS_NOBLE_ABEL_HPP_

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

class NobleAbel : public EosBase<NobleAbel> {
 public:
  NobleAbel() = default;
  PORTABLE_INLINE_FUNCTION
  NobleAbel(Real gm1, Real Cv, Real bb, Real qq,
            const MeanAtomicProperties &AZbar = MeanAtomicProperties())
      : _Cv(Cv), _gm1(gm1), _bb(bb), _qq(qq), _T0(ROOM_TEMPERATURE),
        _P0(ATMOSPHERIC_PRESSURE), _qp(0.0),
        _rho0(robust::ratio(_P0, gm1 * Cv * _T0 + bb * _P0)),
        _vol0(robust::ratio(gm1 * Cv * _T0 + bb * _P0, _P0)), _sie0(Cv * _T0 + qq),
        _bmod0(robust::ratio(_rho0 * Cv * _T0 * gm1 * (gm1 + 1.0),
                             (1.0 - bb * _rho0) * (1.0 - bb * _rho0))),
        _dpde0(robust::ratio(_rho0 * gm1, 1.0 - bb * _rho0)), _AZbar(AZbar) {
    CheckParams();
  }
  PORTABLE_INLINE_FUNCTION
  NobleAbel(Real gm1, Real Cv, Real bb, Real qq, Real qp, Real T0, Real P0,
            const MeanAtomicProperties &AZbar = MeanAtomicProperties())
      : _Cv(Cv), _gm1(gm1), _bb(bb), _qq(qq), _T0(T0), _P0(P0), _qp(qp),
        _rho0(robust::ratio(P0, gm1 * Cv * T0 + bb * P0)),
        _vol0(robust::ratio(gm1 * Cv * T0 + bb * P0, P0)), _sie0(Cv * T0 + qq),
        _bmod0(robust::ratio(_rho0 * Cv * T0 * gm1 * (gm1 + 1.0),
                             (1.0 - bb * _rho0) * (1.0 - bb * _rho0))),
        _dpde0(robust::ratio(_rho0 * gm1, 1.0 - bb * _rho0)), _AZbar(AZbar) {
    CheckParams();
  }
  NobleAbel GetOnDevice() { return *this; }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return std::max(robust::SMALL(), (sie - _qq) / _Cv);
  }
  PORTABLE_INLINE_FUNCTION void CheckParams() const {
    PORTABLE_ALWAYS_REQUIRE(_Cv > 0, "Heat capacity must be positive");
    PORTABLE_ALWAYS_REQUIRE(_gm1 >= 0, "Gruneisen parameter must be positive");
    PORTABLE_ALWAYS_REQUIRE(_bb >= 0, "Covolume must be positive");
    _AZbar.CheckParams();
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return std::max(_qq, _Cv * temperature + _qq);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return std::max(robust::SMALL(),
                    robust::ratio(_gm1 * rho * _Cv * temperature, 1.0 - _bb * rho));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return std::max(robust::SMALL(),
                    robust::ratio(_gm1 * rho * (sie - _qq), 1.0 - _bb * rho));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real MinInternalEnergyFromDensity(
      const Real rho, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return _qq;
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real
  EntropyFromDensityTemperature(const Real rho, const Real temperature,
                                Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real vol = robust::ratio(1.0, rho);
    return _Cv * std::log(robust::ratio(temperature, _T0) + robust::SMALL()) +
           _gm1 * _Cv *
               std::log(robust::ratio(vol - _bb, _vol0 - _bb) + robust::SMALL()) +
           _qp;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real EntropyFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real vol = robust::ratio(1.0, rho);
    return _Cv * std::log(robust::ratio(sie - _qq, _sie0 - _qq) + robust::SMALL()) +
           _gm1 * _Cv *
               std::log(robust::ratio(vol - _bb, _vol0 - _bb) + robust::SMALL()) +
           _qp;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return _Cv;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return _Cv;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return std::max(robust::SMALL(),
                    robust::ratio(_gm1 * (_gm1 + 1.0) * rho * _Cv * temperature,
                                  (1.0 - _bb * rho) * (1.0 - _bb * rho)));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return std::max(robust::SMALL(),
                    robust::ratio(_gm1 * (_gm1 + 1.0) * rho * (sie - _qq),
                                  (1.0 - _bb * rho) * (1.0 - _bb * rho)));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return robust::ratio(_gm1, (1.0 - _bb * rho));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return robust::ratio(_gm1, (1.0 - _bb * rho));
  }
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
    temp = _T0;
    sie = _sie0;
    press = _P0;
    cv = _Cv;
    bmod = _bmod0;
    dpde = _dpde0;
  }
  // Generic functions provided by the base class. These contain e.g. the vector
  // overloads that use the scalar versions declared here
  SG_ADD_DEFAULT_MEAN_ATOMIC_FUNCTIONS(_AZbar)
  SG_ADD_BASE_CLASS_USINGS(NobleAbel)
  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return 0; }
  static constexpr unsigned long PreferredInput() { return _preferred_input; }
  static inline unsigned long scratch_size(std::string method, unsigned int nelements) {
    return 0;
  }
  static inline unsigned long max_scratch_size(unsigned int nelements) { return 0; }
  PORTABLE_INLINE_FUNCTION void PrintParams() const {
    printf("Noble-Abel Parameters:\nGamma = %g\nCv    = %g\nb     = %g\nq     = "
           "%g\n",
           _gm1 + 1.0, _Cv, _bb, _qq);
    _AZbar.PrintParams();
  }
  template <typename Indexer_t>
  PORTABLE_INLINE_FUNCTION void
  DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                       Indexer_t &&lambda, Real &rho, Real &sie) const {
    sie = std::max(_qq, _Cv * temp + _qq);
    rho =
        std::max(robust::SMALL(), robust::ratio(press, _gm1 * _Cv * temp + _bb * press));
  }
  inline void Finalize() {}
  static std::string EosType() { return std::string("NobleAbel"); }
  static std::string EosPyType() { return EosType(); }

 private:
  Real _Cv, _gm1, _bb, _qq;
  // optional reference state variables
  Real _T0, _P0, _qp;
  // reference values
  Real _rho0, _vol0, _sie0, _bmod0, _dpde0;
  MeanAtomicProperties _AZbar;
  // static constexpr const Real _T0 = ROOM_TEMPERATURE;
  // static constexpr const Real _P0 = ATMOSPHERIC_PRESSURE;
  static constexpr const unsigned long _preferred_input =
      thermalqs::density | thermalqs::specific_internal_energy;
};

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION void
NobleAbel::FillEos(Real &rho, Real &temp, Real &sie, Real &press, Real &cv, Real &bmod,
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
    sie = press * (1.0 - _bb * rho) / _gm1 + _qq;
  }
  if (output & thermalqs::pressure) press = PressureFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::temperature)
    temp = TemperatureFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::bulk_modulus)
    bmod = BulkModulusFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::specific_heat)
    cv = SpecificHeatFromDensityInternalEnergy(rho, sie);
}

#ifdef SINGULARITY_INSTANTIATE_CLASSES
SG_ADD_TEMPLATE_EXTERNS(NobleAbel, Real *)
#endif // SINGULARITY_INSTANTIATE_CLASSES

} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_EOS_NOBLE_ABEL_HPP_

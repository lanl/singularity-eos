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

#ifndef _SINGULARITY_EOS_EOS_EOS_CARNAHAN_STARLING_HPP_
#define _SINGULARITY_EOS_EOS_EOS_CARNAHAN_STARLING_HPP_

// stdlib
#include <cmath>
#include <cstdio>
#include <string>

// Ports-of-call
#include <ports-of-call/portability.hpp>

// Base stuff
#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/eos_error.hpp>
#include <singularity-eos/base/math_utils.hpp>
#include <singularity-eos/base/robust_utils.hpp>
#include <singularity-eos/base/root-finding-1d/root_finding.hpp>
#include <singularity-eos/eos/eos_base.hpp>

namespace singularity {

using namespace eos_base;

class CarnahanStarling : public EosBase<CarnahanStarling> {
 public:
  CarnahanStarling() = default;
  PORTABLE_INLINE_FUNCTION CarnahanStarling(Real gm1, Real Cv, Real bb, Real qq)
      : _Cv(Cv), _gm1(gm1), _bb(bb), _qq(qq), _T0(ROOM_TEMPERATURE),
        _P0(ATMOSPHERIC_PRESSURE), _qp(0.0),
        _rho0(DensityFromPressureTemperature(_P0, _T0)), _vol0(robust::ratio(1.0, _rho0)),
        _sie0(Cv * _T0 + qq),
        _bmod0(_rho0 * Cv * _T0 *
               (PartialRhoZedFromDensity(_rho0) +
                ZedFromDensity(_rho0) * ZedFromDensity(_rho0) * gm1)),
        _dpde0(_rho0 * ZedFromDensity(_rho0) * gm1) {
    checkParams();
  }
  PORTABLE_INLINE_FUNCTION CarnahanStarling(Real gm1, Real Cv, Real bb, Real qq, Real qp,
                                            Real T0, Real P0)
      : _Cv(Cv), _gm1(gm1), _bb(bb), _qq(qq), _T0(T0), _P0(P0), _qp(qp),
        _rho0(DensityFromPressureTemperature(P0, T0)), _vol0(robust::ratio(1.0, _rho0)),
        _sie0(Cv * T0 + qq),
        _bmod0(_rho0 * Cv * T0 *
               (PartialRhoZedFromDensity(_rho0) +
                ZedFromDensity(_rho0) * ZedFromDensity(_rho0) * gm1)),
        _dpde0(_rho0 * ZedFromDensity(_rho0) * gm1) {
    checkParams();
  }
  CarnahanStarling GetOnDevice() { return *this; }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return std::max(robust::SMALL(), (sie - _qq) / _Cv);
  }
  PORTABLE_INLINE_FUNCTION void checkParams() const {
    PORTABLE_ALWAYS_REQUIRE(_Cv >= 0, "Heat capacity must be positive");
    PORTABLE_ALWAYS_REQUIRE(_gm1 >= 0, "Gruneisen parameter must be positive");
    PORTABLE_ALWAYS_REQUIRE(_bb >= 0, "Covolume must be positive");
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real ZedFromDensity(
      const Real rho, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real eta = _bb * rho;
    const Real zed =
        1.0 + robust::ratio(eta * (4.0 - 2.0 * eta), math_utils::pow<3>(1.0 - eta));
    return zed;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PartialRhoZedFromDensity(
      const Real rho, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real eta = _bb * rho;
    return 1.0 + robust::ratio(eta * (8.0 - 2.0 * eta), math_utils::pow<4>(1.0 - eta));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real DensityFromPressureTemperature(
      const Real press, const Real temperature, const Real guess = robust::SMALL(),
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    Real real_root;
    auto poly = [=](Real dens) {
      return _Cv * temperature * _gm1 * dens * ZedFromDensity(dens);
    };
    using RootFinding1D::findRoot;
    using RootFinding1D::Status;
    static constexpr Real xtol = 1.0e-12;
    static constexpr Real ytol = 1.0e-12;
    static constexpr Real lo_bound = robust::SMALL();
    const Real hi_bound = robust::ratio(1.0, _bb);
    auto status = findRoot(poly, press, guess, lo_bound, hi_bound, xtol, ytol, real_root);
    if (status != Status::SUCCESS) {
      // Root finder failed even though the solution was bracketed... this is an error
      EOS_ERROR("*** (Warning) DensityFromPressureTemperature :: Convergence not met in "
                "Carnahan-Starling EoS (root finder util) ***\n");
      real_root = -1.0;
    }
    return std::max(robust::SMALL(), real_root);
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
    const Real zed = ZedFromDensity(rho);
    return std::max(robust::SMALL(), zed * rho * temperature * _Cv * _gm1);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real zed = ZedFromDensity(rho);
    return std::max(robust::SMALL(), zed * rho * (sie - _qq) * _gm1);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real
  EntropyFromDensityTemperature(const Real rho, const Real temperature,
                                Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real vol = robust::ratio(1.0, rho);
    const Real one_by_vb = robust::ratio(1.0, vol - _bb);
    const Real one_by_v0b = robust::ratio(1.0, _vol0 - _bb);
    return _Cv * std::log(robust::ratio(temperature, _T0) + robust::SMALL()) +
           _gm1 * _Cv * std::log(robust::ratio(vol, _vol0) + robust::SMALL()) -
           _gm1 * _Cv * _bb * (one_by_vb - one_by_v0b) *
               (4.0 + _bb * (one_by_vb + one_by_v0b)) +
           _qp;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real EntropyFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real vol = robust::ratio(1.0, rho);
    const Real one_by_vb = robust::ratio(1.0, vol - _bb);
    const Real one_by_v0b = robust::ratio(1.0, _vol0 - _bb);
    return _Cv * std::log(robust::ratio(sie - _qq, _sie0 - _qq) + robust::SMALL()) +
           _gm1 * _Cv * std::log(robust::ratio(vol, _vol0) + robust::SMALL()) -
           _gm1 * _Cv * _bb * (one_by_vb - one_by_v0b) *
               (4.0 + _bb * (one_by_vb + one_by_v0b)) +
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
                    rho * _Cv * temperature * _gm1 *
                        (PartialRhoZedFromDensity(rho) +
                         ZedFromDensity(rho) * ZedFromDensity(rho) * _gm1));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return std::max(robust::SMALL(),
                    rho * (sie - _qq) * _gm1 *
                        (PartialRhoZedFromDensity(rho) +
                         ZedFromDensity(rho) * ZedFromDensity(rho) * _gm1));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return ZedFromDensity(rho) * _gm1;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return ZedFromDensity(rho) * _gm1;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real MinInternalEnergyFromDensity(
      const Real rho, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return _qq;
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
  SG_ADD_BASE_CLASS_USINGS(CarnahanStarling)
  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return 0; }
  static constexpr unsigned long PreferredInput() { return _preferred_input; }
  static inline unsigned long scratch_size(std::string method, unsigned int nelements) {
    return 0;
  }
  static inline unsigned long max_scratch_size(unsigned int nelements) { return 0; }
  PORTABLE_INLINE_FUNCTION void PrintParams() const {
    printf("Carnahan-Starling Parameters:\nGamma = %g\nCv    = %g\nb     = %g\nq     = "
           "%g\n",
           _gm1 + 1.0, _Cv, _bb, _qq);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                       Indexer_t &&lambda, Real &rho, Real &sie) const {
    sie = std::max(_qq, _Cv * temp + _qq);
    rho = DensityFromPressureTemperature(press, temp);
  }
  inline void Finalize() {}
  static std::string EosType() { return std::string("CarnahanStarling"); }
  static std::string EosPyType() { return EosType(); }

 private:
  Real _Cv, _gm1, _bb, _qq;
  // optional reference state variables
  Real _T0, _P0, _qp;
  // reference values
  Real _rho0, _vol0, _sie0, _bmod0, _dpde0;
  // static constexpr const Real _T0 = ROOM_TEMPERATURE;
  // static constexpr const Real _P0 = ATMOSPHERIC_PRESSURE;
  static constexpr const unsigned long _preferred_input =
      thermalqs::density | thermalqs::specific_internal_energy;
};
template <typename Indexer_t = Real *>
PORTABLE_INLINE_FUNCTION void CarnahanStarling::FillEos(Real &rho, Real &temp, Real &sie,
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
    sie = robust::ratio(press, ZedFromDensity(rho) * rho * _gm1) + _qq;
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

#endif // _SINGULARITY_EOS_EOS_EOS_CARNAHAN_STARLING_HPP_

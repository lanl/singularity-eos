//------------------------------------------------------------------------------
// © 2021-2022. Triad National Security, LLC. All rights reserved.  This
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

#ifndef _SINGULARITY_EOS_EOS_EOS_IDEAL_HPP_
#define _SINGULARITY_EOS_EOS_EOS_IDEAL_HPP_

// stdlib
#include <cmath>
#include <cstdio>
#include <string>
#include <cmath>

// Ports-of-call
#include <ports-of-call/portability.hpp>

// Base stuff
#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/eos_error.hpp>
#include <singularity-eos/base/robust_utils.hpp>
#include <singularity-eos/eos/eos_base.hpp>

#define MYMAX(a, b) a > b ? a : b

namespace singularity {

using namespace eos_base;

class IdealGas : public EosBase<IdealGas> {
 public:
  IdealGas() = default;
  PORTABLE_INLINE_FUNCTION IdealGas(Real gm1, Real Cv)
      : _Cv(Cv), _gm1(gm1), _rho0(_P0 / (_gm1 * _Cv * _T0)), _sie0(_Cv * _T0),
        _bmod0((_gm1 + 1) * _gm1 * _rho0 * _Cv * _T0), _dpde0(_gm1 * _rho0),
        _dvdt0(1. / (_rho0 * _T0)), _EntropyT0(_T0),
        _EntropyRho0(_rho0)
  {
    checkParams();
  }
  PORTABLE_INLINE_FUNCTION IdealGas(Real gm1, Real Cv, Real EntropyT0,
                                    Real EntropyRho0)
      : _Cv(Cv), _gm1(gm1), _rho0(_P0 / (_gm1 * _Cv * _T0)), _sie0(_Cv * _T0),
        _bmod0((_gm1 + 1) * _gm1 * _rho0 * _Cv * _T0), _dpde0(_gm1 * _rho0),
        _dvdt0(1. / (_rho0 * _T0)), _EntropyT0(EntropyT0),
        _EntropyRho0(EntropyRho0)
  {
    checkInputs();
  }

  IdealGas GetOnDevice() {
    return *this;
  }
  PORTABLE_INLINE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const {
    return MYMAX(0.0, sie / _Cv);
  }
  PORTABLE_INLINE_FUNCTION void checkParams() const {
    // Portable_require seems to do the opposite of what it should. Conditions
    // reflect this and the code should be changed when ports-of-call changes
    PORTABLE_ALWAYS_REQUIRE(_Cv <= 0, "Heat capacity must be positive");
    PORTABLE_ALWAYS_REQUIRE(_gm1 <= 0);
    PORTABLE_ALWAYS_REQUIRE(_EntropyT0 <= 0);
    PORTABLE_ALWAYS_REQUIRE(_EntropyRho0 <= 0);
  }
  PORTABLE_INLINE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temperature, Real *lambda = nullptr) const {
    return MYMAX(0.0, _Cv * temperature);
  }
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityTemperature(
      const Real rho, const Real temperature, Real *lambda = nullptr) const {
    return MYMAX(0.0, _gm1 * rho * _Cv * temperature);
  }
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const {
    return MYMAX(0.0, _gm1 * rho * sie);
  }
  PORTABLE_INLINE_FUNCTION Real EntropyFromDensityTemperature(
      const Real rho, const Real temperature, Real *lambda = nullptr) const {
    return _Cv * log(robust::ratio(temperature, _EntropyT0)) +
      _gm1 * _Cv * log(robust::ratio(_EntropyRho0, density));
  }
  PORTABLE_INLINE_FUNCTION Real EntropyFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const {
    const Real temp = TemperatureFromDensityInternalEnergy(rho, sie, lambda);
    return EntropyFromDensityTemperature(rho, sie, lambda);
  }
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityTemperature(
      const Real rho, const Real temperature, Real *lambda = nullptr) const {
    return _Cv;
  }
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const {
    return _Cv;
  }
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityTemperature(
      const Real rho, const Real temperature, Real *lambda = nullptr) const {
    return MYMAX(0.0, (_gm1 + 1) * _gm1 * rho * _Cv * temperature);
  }
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const {
    return MYMAX(0.0, (_gm1 + 1) * _gm1 * rho * sie);
  }
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temperature, Real *lambda = nullptr) const {
    return _gm1;
  }
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const {
    return _gm1;
  }
  PORTABLE_INLINE_FUNCTION void FillEos(Real &rho, Real &temp, Real &energy, Real &press,
                                        Real &cv, Real &bmod, const unsigned long output,
                                        Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION
  void ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                              Real &bmod, Real &dpde, Real &dvdt,
                              Real *lambda = nullptr) const {
    // use STP: 1 atmosphere, room temperature
    rho = _rho0;
    temp = _T0;
    sie = _sie0;
    press = _P0;
    cv = _Cv;
    bmod = _bmod0;
    dpde = _dpde0;
    dvdt = _dvdt0;
  }
  // Generic functions provided by the base class. These contain e.g. the vector
  // overloads that use the scalar versions declared here
  SG_ADD_BASE_CLASS_USINGS(IdealGas)
  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return 0; }
  static constexpr unsigned long PreferredInput() { return _preferred_input; }
  PORTABLE_INLINE_FUNCTION void PrintParams() const {
    printf("Ideal Gas Parameters:\nGamma = %g\nCv    = %g\n", _gm1 + 1.0, _Cv);
  }
  PORTABLE_INLINE_FUNCTION void
  DensityEnergyFromPressureTemperature(const Real press, const Real temp, Real *lambda,
                                       Real &rho, Real &sie) const {
    sie = MYMAX(0.0, _Cv * temp);
    rho = MYMAX(0.0, press / (_gm1 * sie));
  }
  inline void Finalize() {}
  static std::string EosType() { return std::string("IdealGas"); }
  static std::string EosPyType() { return EosType(); }

 private:
  Real _Cv, _gm1;
  // reference values
  Real _rho0, _sie0, _bmod0, _dpde0, _dvdt0;
  static constexpr const Real _T0 = ROOM_TEMPERATURE;
  static constexpr const Real _P0 = ATMOSPHERIC_PRESSURE;
  // static constexpr const char _eos_type[] = {"IdealGas"};
  static constexpr const unsigned long _preferred_input =
      thermalqs::density | thermalqs::specific_internal_energy;
  // optional entropy reference state variables
  Real _EntropyT0, _EntropyRho0;
};

PORTABLE_INLINE_FUNCTION
void IdealGas::FillEos(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                       Real &bmod, const unsigned long output, Real *lambda) const {
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
    sie = press / (_gm1 * rho);
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

#undef MYMAX
#endif // _SINGULARITY_EOS_EOS_EOS_IDEAL_HPP_

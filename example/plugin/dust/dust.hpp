//------------------------------------------------------------------------------
// © 2024. Triad National Security, LLC. All rights reserved.  This
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

#ifndef _EXAMPLE_PLUGIN_DUST_
#define _EXAMPLE_PLUGIN_DUST_

// Ports-of-call
#include <ports-of-call/portability.hpp>

// Base stuff
#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/error_utils.hpp>
#include <singularity-eos/base/robust_utils.hpp>
#include <singularity-eos/eos/eos_base.hpp>

namespace singularity {

using namespace eos_base;

class Dust : public EosBase<Dust> {
 public:
  PORTABLE_INLINE_FUNCTION Dust() : _Cv(1) {}
  PORTABLE_INLINE_FUNCTION Dust(Real Cv) : _Cv(Cv) {}

  Dust GetOnDevice() { return *this; }
  inline void Finalize() {}

  PORTABLE_INLINE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const {
    return sie / _Cv;
  }
  PORTABLE_INLINE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temperature, Real *lambda = nullptr) {
    return _Cv * temperature;
  }
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityTemperature(
      const Real rho, const Real temperature, Real *lambda = nullptr) const {
    return 0;
  }
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const {
    return 0;
  }
  PORTABLE_INLINE_FUNCTION Real
  MinInternalEnergyFromDensity(const Real rho, Real *lambda = nullptr) const {
    return 0.0;
  };

  PORTABLE_INLINE_FUNCTION Real EntropyFromDensityTemperature(
      const Real rho, const Real temperature, Real *lambda = nullptr) const {
    return 0;
  }
  PORTABLE_INLINE_FUNCTION Real EntropyFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const {
    const Real temp = TemperatureFromDensityInternalEnergy(rho, sie, lambda);
    return 0;
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
    return 0;
  }
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const {
    return 0;
  }
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temperature, Real *lambda = nullptr) const {
    return 0;
  }
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const {
    return 0;
  }
  PORTABLE_INLINE_FUNCTION void FillEos(Real &rho, Real &temp, Real &energy, Real &press,
                                        Real &cv, Real &bmod, const unsigned long output,
                                        Real *lambda = nullptr) const {
    press = 0;
    bmod = 0;
    cv = _Cv;
    if (output & thermalqs::density) {
      UNDEFINED_ERROR;
    }
    if (output & thermalqs::specific_internal_energy) {
      if (!(output & thermalqs::temperature)) {
        UNDEFINED_ERROR;
      }
      energy = _Cv * temp;
    }
    if (output & thermalqs::temperature) {
      if (!(output & thermalqs::specific_internal_energy)) {
        temp = energy / _Cv;
      }
    }
  };
  PORTABLE_INLINE_FUNCTION
  void ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                              Real &bmod, Real &dpde, Real &dvdt,
                              Real *lambda = nullptr) const {
    rho = temp = sie = 1;
    press = bmod = dpde = 0;
    dvdt = 0;
  }
  SG_ADD_BASE_CLASS_USINGS(Dust)

  PORTABLE_INLINE_FUNCTION int nlambda() const noexcept { return 0; }
  static constexpr unsigned long PreferredInput() { return _preferred_input; }
  static inline unsigned long scratch_size(std::string method, unsigned int nelements) {
    return 0;
  }
  static inline unsigned long max_scratch_size(unsigned int nelements) { return 0; }
  PORTABLE_INLINE_FUNCTION void PrintParams() const {
    printf("Dust Parameters:\nCv    = %g\n", _Cv);
  }
  static std::string EosType() { return std::string("Dust"); }
  static std::string EosPyType() { return EosType(); }

 private:
  Real _Cv = 1;
  static constexpr const unsigned long _preferred_input =
      thermalqs::density | thermalqs::specific_internal_energy | thermalqs::temperature;
};

} // namespace singularity

#endif // _EXAMPLE_PLUGIN_DUST_

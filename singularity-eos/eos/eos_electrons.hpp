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

#ifndef _SINGULARITY_EOS_EOS_ELECTRONS
#define _SINGULARITY_EOS_EOS_ELECTRONS

#include <cmath>

// Ports-of-call
#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_errors.hpp>

// Base stuff
#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/robust_utils.hpp>
#include <singularity-eos/eos/eos_base.hpp>

namespace singularity {
using namespace eos_base;

/*
  An ideal gas class designed specifically for electrons in a 3T
  framework. The overall structure of the EOS is still identical to an
  ideal gas:

  P = (gamma - 1) rho_i sie
  sie = Cv T

  but with some caviats.
  1. gamma = 5./3. as for a gas with only 3 translational degrees of freedom
  2. rho_i is the ION mass density. Electron mass density is assumed
     to be negligible.
  3. As such,

  Cv = kb <Z_i> / ((gamma -1) <A> mp)

  where <Z_i> is the average ionization state of the gas and <A> is the
  mean atomic mass of the underlying atomic nuclei.

  <A> is read in at initialization. <Z_i> is passed in via the lambda.
 */
class IdealElectrons : public EosBase<IdealElectrons> {
 public:
  struct Lambda { // Lambda index is ionization state Zi
    enum Index { Zi = 0 };
  };
  static constexpr const Real kb = 1.3806505e-16;  // Boltzmann's constant. In ergs/K
  static constexpr const Real mp = 1.67262171e-24; // proton mass, in g.

  IdealElectrons() = default;
  IdealElectrons(const MeanAtomicProperties &AZbar)
      : _AZbar(AZbar), _Cvbase(robust::ratio(kb, _gm1 * _AZbar.Abar * mp)) {
    CheckParams();
  }

  IdealElectrons GetOnDevice() { return *this; }

  PORTABLE_INLINE_FUNCTION void CheckParams() const {
    _AZbar.CheckParams();
    return;
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return std::max(0.0, sie / _Cv(lambda));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return std::max(0.0, _Cv(lambda) * temperature);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return std::max(0.0, _gm1 * rho * _Cv(lambda) * temperature);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return std::max(0.0, _gm1 * rho * sie);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real MinInternalEnergyFromDensity(
      const Real rho, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return 0.0;
  };

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real
  EntropyFromDensityTemperature(const Real rho, const Real temperature,
                                Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    // TODO(JMM): Is this good enough?
    const Real cv = _Cv(lambda);
    return cv * std::log(robust::ratio(temperature, _T0)) +
           _gm1 * cv * std::log(robust::ratio(_rho0, rho));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real EntropyFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real temp = TemperatureFromDensityInternalEnergy(rho, sie, lambda);
    return EntropyFromDensityTemperature(rho, temp, lambda);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return _Cv(lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return _Cv(lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return std::max(0.0, (_gm1 + 1) * _gm1 * rho * _Cv(lambda) * temperature);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return std::max(0.0, (_gm1 + 1) * _gm1 * rho * sie);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return _gm1;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return _gm1;
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                       Indexer_t &&lambda, Real &rho, Real &sie) const {
    sie = std::max(0.0, _Cv(lambda) * temp);
    rho = std::max(0.0, press / (_gm1 * sie));
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                         Real &bmod, Real &dpde, Real &dvdt,
                         Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    // use STP: 1 atmosphere, room temperature
    const Real Cv = _Cv(lambda);
    rho = _rho0;
    temp = _T0;
    sie = Cv * temp;
    press = _P0;
    cv = Cv;
    bmod = BulkModulusFromDensityTemperature(rho, temp, lambda);
    dpde = _dpde0;
    dvdt = _dvdt0;
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  FillEos(Real &rho, Real &temp, Real &energy, Real &press, Real &cv, Real &bmod,
          const unsigned long output,
          Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;

  // Generic functions provided by the base class. These contain
  // e.g. the vector overloads that use the scalar versions declared
  // here
  SG_ADD_DEFAULT_MEAN_ATOMIC_FUNCTIONS(_AZbar)
  SG_ADD_BASE_CLASS_USINGS(IdealElectrons)
  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return 1; }
  static constexpr unsigned long PreferredInput() { return _preferred_input; }
  static inline unsigned long scratch_size(std::string method, unsigned int nelements) {
    return 0;
  }
  static inline unsigned long max_scratch_size(unsigned int nelements) { return 0; }
  PORTABLE_INLINE_FUNCTION void PrintParams() const {
    printf("Ideal Electrons Parameters:\n");
    _AZbar.PrintParams();
  }
  inline void Finalize() {}
  static std::string EosType() { return std::string("IdealElectrons"); }
  static std::string EosPyType() { return EosType(); }

 private:
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real _Cv(Indexer_t &&lambda) const {
    const Real Z = std::max(lambda[Lambda::Zi], static_cast<Real>(0.0));
    return _Cvbase * Z;
  }

  // TODO(JMM): Change gamma if needed
  // JMM: For some reason some cuda/HIP implementations don't like
  // static constexpr vars. They are not getting properly translated to
  // device, even with --expt-relaxed-constexpr
  // And setting these to const deletes the operator=
  Real _gm1 = (5. / 3.) - 1; // 3 translational DOF
  Real _rho0 = 1.0;
  Real _T0 = ROOM_TEMPERATURE;
  Real _P0 = ATMOSPHERIC_PRESSURE;
  Real _dpde0 = _gm1 * _rho0;
  Real _dvdt0 = 1. / (_rho0 * _T0);
  static constexpr const unsigned long _preferred_input =
      thermalqs::density | thermalqs::specific_internal_energy;
  MeanAtomicProperties _AZbar;
  Real _Cvbase; // cv / <Zi> where Zi is ionization state
};

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION void IdealElectrons::FillEos(Real &rho, Real &temp, Real &sie,
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
    sie = press / (_gm1 * rho);
  }
  if (output & thermalqs::pressure)
    press = PressureFromDensityInternalEnergy(rho, sie, lambda);
  if (output & thermalqs::temperature)
    temp = TemperatureFromDensityInternalEnergy(rho, sie, lambda);
  if (output & thermalqs::bulk_modulus)
    bmod = BulkModulusFromDensityInternalEnergy(rho, sie, lambda);
  if (output & thermalqs::specific_heat)
    cv = SpecificHeatFromDensityInternalEnergy(rho, sie, lambda);
}

} // namespace singularity

#endif //  _SINGULARITY_EOS_EOS_ELECTRONS

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

#ifndef _SINGULARITY_EOS_EOS_EOS_UNITSYSTEM_HPP_
#define _SINGULARITY_EOS_EOS_EOS_UNITSYSTEM_HPP_

#include "stdio.h"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <utility>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_errors.hpp>
#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/eos_error.hpp>
#include <singularity-eos/eos/eos_base.hpp>

namespace singularity {

using namespace eos_base;

// tag dispatch for constructors for UnitSystem
namespace eos_units_init {
[[maybe_unused]] static struct ThermalUnitsInit {
} thermal_units_init_tag;
[[maybe_unused]] static struct LengthTimeUnitsInit {
} length_time_units_init_tag;
} // namespace eos_units_init

template <typename T>
class UnitSystem : public EosBase<UnitSystem<T>> {
 public:
  SG_ADD_BASE_CLASS_USINGS(UnitSystem<T>);
  using BaseType = T;

  // give me std::format or fmt::format...
  static std::string EosType() {
    return std::string("UnitSystem<") + T::EosType() + std::string(">");
  }

  static std::string EosPyType() { return std::string("UnitSystem") + T::EosPyType(); }

  // move semantics ensures dynamic memory comes along for the ride
  // TODO(JMM): Entropy unit needed?
  UnitSystem(T &&t, eos_units_init::ThermalUnitsInit, const Real rho_unit,
             const Real sie_unit, const Real temp_unit)
      : t_(std::forward<T>(t)), rho_unit_(rho_unit), sie_unit_(sie_unit),
        temp_unit_(temp_unit), press_unit_(rho_unit * sie_unit),
        entropy_unit_(sie_unit / temp_unit),
        inv_rho_unit_(1 / rho_unit) // inverses computed to avoid division at runtime
        ,
        inv_sie_unit_(1 / sie_unit), inv_temp_unit_(1 / temp_unit),
        inv_press_unit_(1 / press_unit_), inv_entropy_unit_(1 / entropy_unit_),
        inv_dpde_unit_(sie_unit / press_unit_) // thermo derivatives computed consistently
        ,
        inv_dvdt_unit_(rho_unit * temp_unit) // TODO(JMM): Is this convention weird?
        ,
        inv_dpdr_unit_(rho_unit / press_unit_), inv_dtdr_unit_(rho_unit / temp_unit),
        inv_dtde_unit_(sie_unit / temp_unit) // obviously this is also Cv
        ,
        inv_cv_unit_(temp_unit / sie_unit), inv_bmod_unit_(1 / press_unit_) {
    CheckParams();
  }
  UnitSystem(T &&t, eos_units_init::LengthTimeUnitsInit, const Real time_unit,
             const Real mass_unit, const Real length_unit, const Real temp_unit)
      : UnitSystem(std::forward<T>(t), eos_units_init::thermal_units_init_tag,
                   mass_unit / (length_unit * length_unit * length_unit),
                   length_unit * length_unit / (time_unit * time_unit), temp_unit) {}
  UnitSystem(T &&t, const Real rho_unit, const Real sie_unit, const Real temp_unit)
      : UnitSystem(std::forward<T>(t), eos_units_init::thermal_units_init_tag, rho_unit,
                   sie_unit, temp_unit) {}
  UnitSystem() = default;

  PORTABLE_INLINE_FUNCTION void CheckParams() const {
    PORTABLE_ALWAYS_REQUIRE(rho_unit_ > 0, "Nonzero density unit");
    PORTABLE_ALWAYS_REQUIRE(sie_unit_ > 0, "Nonzero energy unit");
    PORTABLE_ALWAYS_REQUIRE(temp_unit_ > 0, "Nonzero temperature unit");
    t_.CheckParams();
  }

  auto GetOnDevice() {
    return UnitSystem<T>(t_.GetOnDevice(), eos_units_init::thermal_units_init_tag,
                         rho_unit_, sie_unit_, temp_unit_);
  }
  inline void Finalize() { t_.Finalize(); }

  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    const Real temp =
        t_.TemperatureFromDensityInternalEnergy(rho * rho_unit_, sie * sie_unit_, lambda);
    return inv_temp_unit_ * temp;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temperature, Indexer_t &&lambda = nullptr) const {
    const Real sie = t_.InternalEnergyFromDensityTemperature(
        rho * rho_unit_, temperature * temp_unit_, lambda);
    return inv_sie_unit_ * sie;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real PressureFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    const Real P =
        t_.PressureFromDensityInternalEnergy(rho * rho_unit_, sie * sie_unit_, lambda);
    return inv_press_unit_ * P;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real
  MinInternalEnergyFromDensity(const Real rho, Indexer_t &&lambda = nullptr) const {
    const Real S = t_.MinInternalEnergyFromDensity(rho * rho_unit_, lambda);
    return inv_sie_unit_ * S;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real EntropyFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    const Real S =
        t_.EntropyFromDensityInternalEnergy(rho * rho_unit_, sie * sie_unit_, lambda);
    return inv_entropy_unit_ * S;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    const Real cv = t_.SpecificHeatFromDensityInternalEnergy(rho * rho_unit_,
                                                             sie * sie_unit_, lambda);
    return inv_cv_unit_ * cv;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    const Real bmod =
        t_.BulkModulusFromDensityInternalEnergy(rho * rho_unit_, sie * sie_unit_, lambda);
    return inv_bmod_unit_ * bmod;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    const Real gm1 = t_.GruneisenParamFromDensityInternalEnergy(rho * rho_unit_,
                                                                sie * sie_unit_, lambda);
    return gm1;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real PressureFromDensityTemperature(
      const Real rho, const Real temp, Indexer_t &&lambda = nullptr) const {
    const Real P =
        t_.PressureFromDensityTemperature(rho * rho_unit_, temp * temp_unit_, lambda);
    return inv_press_unit_ * P;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real EntropyFromDensityTemperature(
      const Real rho, const Real temp, Indexer_t &&lambda = nullptr) const {
    const Real S =
        t_.EntropyFromDensityTemperature(rho * rho_unit_, temp * temp_unit_, lambda);
    return inv_entropy_unit_ * S;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real SpecificHeatFromDensityTemperature(
      const Real rho, const Real temp, Indexer_t &&lambda = nullptr) const {
    const Real cv =
        t_.SpecificHeatFromDensityTemperature(rho * rho_unit_, temp * temp_unit_, lambda);
    return inv_cv_unit_ * cv;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real BulkModulusFromDensityTemperature(
      const Real rho, const Real temp, Indexer_t &&lambda = nullptr) const {
    const Real bmod =
        t_.BulkModulusFromDensityTemperature(rho * rho_unit_, temp * temp_unit_, lambda);
    return inv_bmod_unit_ * bmod;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temp, Indexer_t &&lambda = nullptr) const {
    const Real gm1 = t_.GruneisenParamFromDensityTemperature(rho * rho_unit_,
                                                             temp * temp_unit_, lambda);
    return gm1;
  }

  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION void
  DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                       Indexer_t &&lambda, Real &rho, Real &sie) const {
    t_.DensityEnergyFromPressureTemperature(press * press_unit_, temp * temp_unit_,
                                            lambda, rho, sie);
    rho *= inv_rho_unit_;
    sie *= inv_sie_unit_;
  }

  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION void FillEos(Real &rho, Real &temp, Real &energy, Real &press,
                                 Real &cv, Real &bmod, const unsigned long output,
                                 Indexer_t &&lambda = nullptr) const {
    // TODO(JMM): Is this general enough? Do I need more switches/scales?
    Real srho = rho_unit_ * rho;
    switch (t_.PreferredInput()) {
    case thermalqs::density | thermalqs::temperature: {
      Real sT = temp_unit_ * temp;
      t_.FillEos(srho, sT, energy, press, cv, bmod, output, lambda);
      energy *= inv_sie_unit_;
      break;
    }
    case thermalqs::density | thermalqs::specific_internal_energy: {
      Real ssie = sie_unit_ * energy;
      t_.FillEos(srho, temp, ssie, press, cv, bmod, output, lambda);
      break;
    }
    default: {
      EOS_ERROR("Didn't find a valid input for ScaledEOS::FillEOS\n");
    }
    }
    press *= inv_press_unit_;
    cv *= inv_cv_unit_;
    bmod *= inv_bmod_unit_;
  }

  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION void ValuesAtReferenceState(Real &rho, Real &temp, Real &sie,
                                                Real &press, Real &cv, Real &bmod,
                                                Real &dpde, Real &dvdt,
                                                Indexer_t &&lambda = nullptr) const {
    t_.ValuesAtReferenceState(rho, temp, sie, press, cv, bmod, dpde, dvdt, lambda);
    rho *= inv_rho_unit_;
    temp *= inv_temp_unit_;
    sie *= inv_sie_unit_;
    press *= inv_press_unit_;
    cv *= inv_cv_unit_;
    bmod *= inv_bmod_unit_;
    dpde *= inv_dpde_unit_;
    dvdt *= inv_dvdt_unit_;
  }

  PORTABLE_FORCEINLINE_FUNCTION Real MinimumDensity() const {
    return inv_rho_unit_ * t_.MinimumDensity();
  }
  PORTABLE_FORCEINLINE_FUNCTION Real MinimumTemperature() const {
    return inv_temp_unit_ * t_.MinimumTemperature();
  }
  PORTABLE_FORCEINLINE_FUNCTION Real MaximumDensity() const {
    return inv_rho_unit_ * t_.MaximumDensity();
  }
  PORTABLE_FORCEINLINE_FUNCTION Real MinimumPressure() const {
    return inv_press_unit_ * t_.MinimumPressure();
  }
  PORTABLE_FORCEINLINE_FUNCTION Real
  MaximumPressureFromTemperature(const Real temp) const {
    return inv_press_unit_ * t_.MaximumPressureFromTemperature(temp_unit_ * temp);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real MeanAtomicMassFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return t_.MeanAtomicMassFromDensityTemperature(rho * rho_unit_,
                                                   temperature * temp_unit_, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real MeanAtomicNumberFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return t_.MeanAtomicNumberFromDensityTemperature(rho * rho_unit_,
                                                     temperature * temp_unit_, lambda);
  }

  // vector implementations
  template <typename LambdaIndexer>
  inline void TemperatureFromDensityInternalEnergy(
      const Real *rhos, const Real *sies, Real *temperatures, Real *scratch,
      const int num, LambdaIndexer &&lambdas, Transform &&transform = Transform()) const {
    transform.x.apply(rho_unit_);
    transform.y.apply(sie_unit_);
    transform.f.apply(inv_temp_unit_);
    t_.TemperatureFromDensityInternalEnergy(rhos, sies, temperatures, scratch, num,
                                            std::forward<LambdaIndexer>(lambdas),
                                            std::forward<Transform>(transform));
  }

  template <typename LambdaIndexer>
  inline void PressureFromDensityTemperature(const Real *rhos, const Real *temperatures,
                                             Real *pressures, Real *scratch,
                                             const int num, LambdaIndexer &&lambdas,
                                             Transform &&transform = Transform()) const {
    transform.x.apply(rho_unit_);
    transform.y.apply(temp_unit_);
    transform.f.apply(inv_press_unit_);
    t_.PressureFromDensityTemperature(rhos, temperatures, pressures, scratch, num,
                                      std::forward<LambdaIndexer>(lambdas),
                                      std::forward<Transform>(transform));
  }

  template <typename LambdaIndexer>
  inline void
  PressureFromDensityInternalEnergy(const Real *rhos, const Real *sies, Real *pressures,
                                    Real *scratch, const int num, LambdaIndexer &&lambdas,
                                    Transform &&transform = Transform()) const {
    transform.x.apply(rho_unit_);
    transform.y.apply(sie_unit_);
    transform.f.apply(inv_press_unit_);
    t_.PressureFromDensityInternalEnergy(rhos, sies, pressures, scratch, num,
                                         std::forward<LambdaIndexer>(lambdas),
                                         std::forward<Transform>(transform));
  }

  template <typename LambdaIndexer>
  inline void MinInternalEnergyFromDensity(const Real *rhos, Real *sies, Real *scratch,
                                           const int num, LambdaIndexer &&lambdas,
                                           Transform &&transform = Transform()) const {
    transform.x.apply(rho_unit_);
    transform.f.apply(sie_unit_);
    t_.MinInternalEnergyFromDensity(rhos, sies, scratch, num,
                                    std::forward<LambdaIndexer>(lambdas),
                                    std::forward<Transform>(transform));
  }

  template <typename LambdaIndexer>
  inline void SpecificHeatFromDensityTemperature(
      const Real *rhos, const Real *temperatures, Real *cvs, Real *scratch, const int num,
      LambdaIndexer &&lambdas, Transform &&transform = Transform()) const {
    transform.x.apply(rho_unit_);
    transform.y.apply(temp_unit_);
    transform.f.apply(inv_cv_unit_);
    t_.SpecificHeatFromDensityTemperature(rhos, temperatures, cvs, scratch, num,
                                          std::forward<LambdaIndexer>(lambdas),
                                          std::forward<Transform>(transform));
  }

  template <typename LambdaIndexer>
  inline void SpecificHeatFromDensityInternalEnergy(
      const Real *rhos, const Real *sies, Real *cvs, Real *scratch, const int num,
      LambdaIndexer &&lambdas, Transform &&transform = Transform()) const {
    transform.x.apply(rho_unit_);
    transform.y.apply(sie_unit_);
    transform.f.apply(inv_cv_unit_);
    t_.SpecificHeatFromDensityInternalEnergy(rhos, sies, cvs, scratch, num,
                                             std::forward<LambdaIndexer>(lambdas),
                                             std::forward<Transform>(transform));
  }

  template <typename LambdaIndexer>
  inline void BulkModulusFromDensityTemperature(
      const Real *rhos, const Real *temperatures, Real *bmods, Real *scratch,
      const int num, LambdaIndexer &&lambdas, Transform &&transform = Transform()) const {
    transform.x.apply(rho_unit_);
    transform.y.apply(temp_unit_);
    transform.f.apply(inv_bmod_unit_);
    t_.BulkModulusFromDensityTemperature(rhos, temperatures, bmods, scratch, num,
                                         std::forward<LambdaIndexer>(lambdas),
                                         std::forward<Transform>(transform));
  }

  template <typename LambdaIndexer>
  inline void BulkModulusFromDensityInternalEnergy(
      const Real *rhos, const Real *sies, Real *bmods, Real *scratch, const int num,
      LambdaIndexer &&lambdas, Transform &&transform = Transform()) const {
    transform.x.apply(rho_unit_);
    transform.y.apply(sie_unit_);
    transform.f.apply(inv_bmod_unit_);
    t_.BulkModulusFromDensityInternalEnergy(rhos, sies, bmods, scratch, num,
                                            std::forward<LambdaIndexer>(lambdas),
                                            std::forward<Transform>(transform));
  }

  template <typename LambdaIndexer>
  inline void GruneisenParamFromDensityTemperature(
      const Real *rhos, const Real *temperatures, Real *gm1s, Real *scratch,
      const int num, LambdaIndexer &&lambdas, Transform &&transform = Transform()) const {
    transform.x.apply(rho_unit_);
    transform.y.apply(temp_unit_);
    t_.GruneisenParamFromDensityTemperature(rhos, temperatures, gm1s, scratch, num,
                                            std::forward<LambdaIndexer>(lambdas),
                                            std::forward<Transform>(transform));
  }

  template <typename LambdaIndexer>
  inline void GruneisenParamFromDensityInternalEnergy(
      const Real *rhos, const Real *sies, Real *gm1s, Real *scratch, const int num,
      LambdaIndexer &&lambdas, Transform &&transform = Transform()) const {
    transform.x.apply(rho_unit_);
    transform.y.apply(sie_unit_);
    t_.GruneisenParamFromDensityInternalEnergy(rhos, sies, gm1s, scratch, num,
                                               std::forward<LambdaIndexer>(lambdas),
                                               std::forward<Transform>(transform));
  }

  template <typename LambdaIndexer>
  inline void InternalEnergyFromDensityTemperature(
      const Real *rhos, const Real *temperatures, Real *sies, Real *scratch,
      const int num, LambdaIndexer &&lambdas, Transform &&transform = Transform()) const {
    transform.x.apply(rho_unit_);
    transform.y.apply(temp_unit_);
    transform.f.apply(inv_sie_unit_);
    t_.InternalEnergyFromDensityTemperature(rhos, temperatures, sies, scratch, num,
                                            std::forward<LambdaIndexer>(lambdas),
                                            std::forward<Transform>(transform));
  }

  template <typename LambdaIndexer>
  inline void EntropyFromDensityTemperature(const Real *rhos, const Real *temperatures,
                                            Real *entropies, Real *scratch, const int num,
                                            LambdaIndexer &&lambdas,
                                            Transform &&transform = Transform()) const {
    transform.x.apply(rho_unit_);
    transform.y.apply(temp_unit_);
    transform.f.apply(inv_entropy_unit_);
    t_.EntropyFromDensityTemperature(rhos, temperatures, entropies, scratch, num,
                                     std::forward<LambdaIndexer>(lambdas),
                                     std::forward<Transform>(transform));
  }

  template <typename LambdaIndexer>
  inline void
  EntropyFromDensityInternalEnergy(const Real *rhos, const Real *sies, Real *entropies,
                                   Real *scratch, const int num, LambdaIndexer &&lambdas,
                                   Transform &&transform = Transform()) const {
    transform.x.apply(rho_unit_);
    transform.y.apply(sie_unit_);
    transform.f.apply(inv_entropy_unit_);
    t_.EntropyFromDensityInternalEnergy(rhos, sies, entropies, scratch, num,
                                        std::forward<LambdaIndexer>(lambdas),
                                        std::forward<Transform>(transform));
  }

  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return t_.nlambda(); }

  static constexpr unsigned long PreferredInput() { return T::PreferredInput(); }

  static inline unsigned long scratch_size(std::string method, unsigned int nelements) {
    return T::scratch_size(method, nelements);
  }

  static inline unsigned long max_scratch_size(unsigned int nelements) {
    return T::max_scratch_size(nelements);
  }

  PORTABLE_FUNCTION void PrintParams() const {
    t_.PrintParams();
    printf("Units = %e %e %e %e\n", rho_unit_, sie_unit_, temp_unit_, press_unit_);
  }

  SG_ADD_MODIFIER_METHODS(T, t_);
  SG_ADD_MODIFIER_MEAN_METHODS(t_)

 private:
  T t_;

  // JMM: The compiler deletes GetOnDevice if I make these const. So... they're not
  Real rho_unit_, sie_unit_, temp_unit_, press_unit_, entropy_unit_;
  Real inv_rho_unit_, inv_sie_unit_, inv_temp_unit_, inv_press_unit_, inv_entropy_unit_;
  Real inv_dpde_unit_, inv_dvdt_unit_, inv_dpdr_unit_, inv_dtdr_unit_, inv_dtde_unit_;
  Real inv_cv_unit_, inv_bmod_unit_;
};

} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_EOS_UNITSYSTEM_HPP_

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

#ifndef _SINGULARITY_EOS_EOS_SHIFTED_EOS_
#define _SINGULARITY_EOS_EOS_SHIFTED_EOS_

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

template <typename T>
class ShiftedEOS : public EosBase<ShiftedEOS<T>> {
 public:
  SG_ADD_BASE_CLASS_USINGS(ShiftedEOS<T>);
  using BaseType = T;

  // give me std::format or fmt::format...
  static std::string EosType() {
    return std::string("ShiftedEOS<") + T::EosType() + std::string(">");
  }

  static std::string EosPyType() { return std::string("Shifted") + T::EosPyType(); }

  // move semantics ensures dynamic memory comes along for the ride
  ShiftedEOS(T &&t, const Real shift) : t_(std::forward<T>(t)), shift_(shift) {
    CheckParams();
  }
  ShiftedEOS() = default;

  PORTABLE_INLINE_FUNCTION void CheckParams() const {
    PORTABLE_ALWAYS_REQUIRE(!std::isnan(shift_), "Shift must be a number");
    t_.CheckParams();
  }

  auto GetOnDevice() { return ShiftedEOS<T>(t_.GetOnDevice(), shift_); }
  inline void Finalize() { t_.Finalize(); }

  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    return t_.TemperatureFromDensityInternalEnergy(rho, sie - shift_, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temperature, Indexer_t &&lambda = nullptr) const {
    Real energy = t_.InternalEnergyFromDensityTemperature(rho, temperature, lambda);
    return energy + shift_;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real PressureFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    return t_.PressureFromDensityInternalEnergy(rho, sie - shift_, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real
  MinInternalEnergyFromDensity(const Real rho, Indexer_t &&lambda = nullptr) const {
    return t_.MinInternalEnergyFromDensity(rho, lambda) + shift_;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real EntropyFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    return t_.EntropyFromDensityInternalEnergy(rho, sie - shift_, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    return t_.SpecificHeatFromDensityInternalEnergy(rho, sie - shift_, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    return t_.BulkModulusFromDensityInternalEnergy(rho, sie - shift_, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    return t_.GruneisenParamFromDensityInternalEnergy(rho, sie - shift_, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real PressureFromDensityTemperature(
      const Real rho, const Real temperature, Indexer_t &&lambda = nullptr) const {
    return t_.PressureFromDensityTemperature(rho, temperature, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real EntropyFromDensityTemperature(
      const Real rho, const Real temperature, Indexer_t &&lambda = nullptr) const {
    return t_.EntropyFromDensityTemperature(rho, temperature, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real SpecificHeatFromDensityTemperature(
      const Real rho, const Real temperature, Indexer_t &&lambda = nullptr) const {
    return t_.SpecificHeatFromDensityTemperature(rho, temperature, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real BulkModulusFromDensityTemperature(
      const Real rho, const Real temperature, Indexer_t &&lambda = nullptr) const {
    return t_.BulkModulusFromDensityTemperature(rho, temperature, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temperature, Indexer_t &&lambda = nullptr) const {
    return t_.GruneisenParamFromDensityTemperature(rho, temperature, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION void FillEos(Real &rho, Real &temp, Real &energy, Real &press,
                                 Real &cv, Real &bmod, const unsigned long output,
                                 Indexer_t &&lambda = nullptr) const {
    Real senergy;
    switch (t_.PreferredInput()) {
    case thermalqs::density | thermalqs::temperature:
      t_.FillEos(rho, temp, energy, press, cv, bmod, output, lambda);
      energy = energy + shift_;
      break;
    case thermalqs::density | thermalqs::specific_internal_energy:
      senergy = energy - shift_;
      t_.FillEos(rho, temp, senergy, press, cv, bmod, output, lambda);
      break;
    default:
      EOS_ERROR("Didn't find a valid input for ShiftedEOS::FillEOS\n");
    }
  }

  // vector implementations
  inline void shift_sies(const Real *sies, Real *shifted, const int num) const {
    static auto const name =
        singularity::mfuncname::member_func_name(typeid(ShiftedEOS<T>).name(), __func__);
    static auto const cname = name.c_str();
    const auto shift_val = shift_;
    portableFor(
        cname, 0, num,
        PORTABLE_LAMBDA(const int i) { shifted[i] = sies[i] - shift_val; });
  }

  inline void unshift_sies(Real *sies, const int num) const {
    static auto const name =
        singularity::mfuncname::member_func_name(typeid(ShiftedEOS<T>).name(), __func__);
    static auto const cname = name.c_str();
    const auto shift_val = shift_;
    portableFor(
        cname, 0, num, PORTABLE_LAMBDA(const int i) { sies[i] += shift_val; });
  }

  template <typename LambdaIndexer>
  inline void TemperatureFromDensityInternalEnergy(
      const Real *rhos, const Real *sies, Real *temperatures, Real *scratch,
      const int num, LambdaIndexer &&lambdas, Transform &&transform = Transform()) const {
    Real *shifted_sies = scratch;
    shift_sies(sies, shifted_sies, num);
    t_.TemperatureFromDensityInternalEnergy(
        rhos, shifted_sies, temperatures, &scratch[num], num,
        std::forward<LambdaIndexer>(lambdas), std::forward<Transform>(transform));
  }

  template <typename LambdaIndexer>
  inline void PressureFromDensityTemperature(const Real *rhos, const Real *temperatures,
                                             Real *pressures, Real *scratch,
                                             const int num, LambdaIndexer &&lambdas,
                                             Transform &&transform = Transform()) const {
    t_.PressureFromDensityTemperature(rhos, temperatures, pressures, scratch, num,
                                      std::forward<LambdaIndexer>(lambdas),
                                      std::forward<Transform>(transform));
  }

  template <typename LambdaIndexer>
  inline void
  PressureFromDensityInternalEnergy(const Real *rhos, const Real *sies, Real *pressures,
                                    Real *scratch, const int num, LambdaIndexer &&lambdas,
                                    Transform &&transform = Transform()) const {
    Real *shifted_sies = scratch;
    shift_sies(sies, shifted_sies, num);
    t_.PressureFromDensityInternalEnergy(rhos, shifted_sies, pressures, &scratch[num],
                                         num, std::forward<LambdaIndexer>(lambdas),
                                         std::forward<Transform>(transform));
  }

  template <typename LambdaIndexer>
  inline void MinInternalEnergyFromDensity(const Real *rhos, Real *sies, Real *scratch,
                                           const int num, LambdaIndexer &&lambdas,
                                           Transform &&transform = Transform()) const {
    t_.MinInternalEnergyFromDensity(rhos, sies, &scratch[num], num,
                                    std::forward<LambdaIndexer>(lambdas),
                                    std::forward<Transform>(transform));
    unshift_sies(sies, num);
  }

  template <typename LambdaIndexer>
  inline void SpecificHeatFromDensityTemperature(
      const Real *rhos, const Real *temperatures, Real *cvs, Real *scratch, const int num,
      LambdaIndexer &&lambdas, Transform &&transform = Transform()) const {
    t_.SpecificHeatFromDensityTemperature(rhos, temperatures, cvs, scratch, num,
                                          std::forward<LambdaIndexer>(lambdas),
                                          std::forward<Transform>(transform));
  }

  template <typename LambdaIndexer>
  inline void SpecificHeatFromDensityInternalEnergy(
      const Real *rhos, const Real *sies, Real *cvs, Real *scratch, const int num,
      LambdaIndexer &&lambdas, Transform &&transform = Transform()) const {
    Real *shifted_sies = scratch;
    shift_sies(sies, shifted_sies, num);
    t_.SpecificHeatFromDensityInternalEnergy(rhos, shifted_sies, cvs, &scratch[num], num,
                                             std::forward<LambdaIndexer>(lambdas),
                                             std::forward<Transform>(transform));
  }

  template <typename LambdaIndexer>
  inline void BulkModulusFromDensityTemperature(
      const Real *rhos, const Real *temperatures, Real *bmods, Real *scratch,
      const int num, LambdaIndexer &&lambdas, Transform &&transform = Transform()) const {
    t_.BulkModulusFromDensityTemperature(rhos, temperatures, bmods, scratch, num,
                                         std::forward<LambdaIndexer>(lambdas),
                                         std::forward<Transform>(transform));
  }

  template <typename LambdaIndexer>
  inline void BulkModulusFromDensityInternalEnergy(
      const Real *rhos, const Real *sies, Real *bmods, Real *scratch, const int num,
      LambdaIndexer &&lambdas, Transform &&transform = Transform()) const {
    Real *shifted_sies = scratch;
    shift_sies(sies, shifted_sies, num);
    t_.BulkModulusFromDensityInternalEnergy(rhos, shifted_sies, bmods, &scratch[num], num,
                                            std::forward<LambdaIndexer>(lambdas),
                                            std::forward<Transform>(transform));
  }

  template <typename LambdaIndexer>
  inline void GruneisenParamFromDensityTemperature(
      const Real *rhos, const Real *temperatures, Real *gm1s, Real *scratch,
      const int num, LambdaIndexer &&lambdas, Transform &&transform = Transform()) const {
    t_.GruneisenParamFromDensityTemperature(rhos, temperatures, gm1s, scratch, num,
                                            std::forward<LambdaIndexer>(lambdas),
                                            std::forward<Transform>(transform));
  }

  template <typename LambdaIndexer>
  inline void GruneisenParamFromDensityInternalEnergy(
      const Real *rhos, const Real *sies, Real *gm1s, Real *scratch, const int num,
      LambdaIndexer &&lambdas, Transform &&transform = Transform()) const {
    Real *shifted_sies = scratch;
    shift_sies(sies, shifted_sies, num);
    t_.GruneisenParamFromDensityInternalEnergy(rhos, shifted_sies, gm1s, &scratch[num],
                                               num, std::forward<LambdaIndexer>(lambdas),
                                               std::forward<Transform>(transform));
  }

  template <typename LambdaIndexer>
  inline void InternalEnergyFromDensityTemperature(
      const Real *rhos, const Real *temperatures, Real *sies, Real *scratch,
      const int num, LambdaIndexer &&lambdas, Transform &&transform = Transform()) const {
    t_.InternalEnergyFromDensityTemperature(rhos, temperatures, sies, scratch, num,
                                            std::forward<LambdaIndexer>(lambdas),
                                            std::forward<Transform>(transform));
    unshift_sies(sies, num);
  }

  template <typename LambdaIndexer>
  inline void EntropyFromDensityTemperature(const Real *rhos, const Real *temperatures,
                                            Real *entropies, Real *scratch, const int num,
                                            LambdaIndexer &&lambdas,
                                            Transform &&transform = Transform()) const {
    t_.EntropyFromDensityTemperature(rhos, temperatures, entropies, scratch, num,
                                     std::forward<LambdaIndexer>(lambdas),
                                     std::forward<Transform>(transform));
  }

  template <typename LambdaIndexer>
  inline void
  EntropyFromDensityInternalEnergy(const Real *rhos, const Real *sies, Real *entropies,
                                   Real *scratch, const int num, LambdaIndexer &&lambdas,
                                   Transform &&transform = Transform()) const {
    Real *shifted_sies = scratch;
    shift_sies(sies, shifted_sies, num);
    t_.EntropyFromDensityInternalEnergy(rhos, sies, entropies, &scratch[num], num,
                                        std::forward<LambdaIndexer>(lambdas),
                                        std::forward<Transform>(transform));
  }

  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return t_.nlambda(); }

  static constexpr unsigned long PreferredInput() { return T::PreferredInput(); }

  static inline unsigned long scratch_size(std::string method, unsigned int nelements) {
    constexpr char suffix[] = "FromDensityInternalEnergy";
    if (method.rfind(suffix) == method.length() - sizeof(suffix) + 1) {
      return T::scratch_size(method, nelements) + (sizeof(Real) * nelements);
    }
    return T::scratch_size(method, nelements);
  }

  static inline unsigned long max_scratch_size(unsigned int nelements) {
    unsigned long m = T::max_scratch_size(nelements);
    for (auto &&meth :
         {"TemperatureFromDensityInternalEnergy", "PressureFromDensityInternalEnergy",
          "SpecificHeatFromDensityInternalEnergy", "BulkModulusFromDensityInternalEnergy",
          "GruneisenParamFromDensityInternalEnergy"}) {
      m = std::max(m, scratch_size(meth, nelements));
    }
    return m;
  }

  PORTABLE_FUNCTION void PrintParams() const {
    t_.PrintParams();
    printf("shift_value = %f\n", shift_);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION void
  DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                       Indexer_t &&lambda, Real &rho, Real &sie) const {
    t_.DensityEnergyFromPressureTemperature(press, temp, lambda, rho, sie);
    sie = sie + shift_;
  }

  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION void ValuesAtReferenceState(Real &rho, Real &temp, Real &sie,
                                                Real &press, Real &cv, Real &bmod,
                                                Real &dpde, Real &dvdt,
                                                Indexer_t &&lambda = nullptr) const {
    t_.ValuesAtReferenceState(rho, temp, sie, press, cv, bmod, dpde, dvdt, lambda);
    sie += shift_;
  }

  PORTABLE_FORCEINLINE_FUNCTION Real MinimumDensity() const {
    return t_.MinimumDensity();
  }
  PORTABLE_FORCEINLINE_FUNCTION Real MinimumTemperature() const {
    return t_.MinimumTemperature();
  }

  SG_ADD_MODIFIER_METHODS(T, t_);

 private:
  T t_;
  double shift_;
};

} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_SHIFTED_EOS_

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

#ifndef _SINGULARITY_EOS_EOS_MINIMUM_ENERGY_
#define _SINGULARITY_EOS_EOS_MINIMUM_ENERGY_

#include "stdio.h"
#include <cstdlib>
#include <iostream>

#include <ports-of-call/portability.hpp>
#include <singularity-eos/base/eos_error.hpp>
#include <singularity-eos/eos/eos_base.hpp>

namespace singularity {

using namespace eos_base;

template <typename T>
class MinimumEnergy : public EosBase<MinimumEnergy<T>> {
 public:
  // Generic functions provided by the base class. These contain
  // e.g. the vector overloads that use the scalar versions declared
  // here We explicitly list, rather than using the macro because we
  // overload some methods.

  // TODO(JMM): The modifier EOS's should probably call the specific
  // sub-functions of the class they modify so that they can leverage,
  // e.g., an especially performant or special version of these
  using EosBase<MinimumEnergy<T>>::TemperatureFromDensityInternalEnergy;
  using EosBase<MinimumEnergy<T>>::InternalEnergyFromDensityTemperature;
  using EosBase<MinimumEnergy<T>>::PressureFromDensityTemperature;
  using EosBase<MinimumEnergy<T>>::PressureFromDensityInternalEnergy;
  using EosBase<MinimumEnergy<T>>::MinInternalEnergyFromDensity;
  using EosBase<MinimumEnergy<T>>::EntropyFromDensityTemperature;
  using EosBase<MinimumEnergy<T>>::EntropyFromDensityInternalEnergy;
  using EosBase<MinimumEnergy<T>>::SpecificHeatFromDensityTemperature;
  using EosBase<MinimumEnergy<T>>::SpecificHeatFromDensityInternalEnergy;
  using EosBase<MinimumEnergy<T>>::BulkModulusFromDensityTemperature;
  using EosBase<MinimumEnergy<T>>::BulkModulusFromDensityInternalEnergy;
  using EosBase<MinimumEnergy<T>>::GruneisenParamFromDensityTemperature;
  using EosBase<MinimumEnergy<T>>::GruneisenParamFromDensityInternalEnergy;
  using EosBase<MinimumEnergy<T>>::FillEos;

  using BaseType = T;

  // give me std::format or fmt::format...
  static std::string EosType() {
    return std::string("MinimumEnergy<") + T::EosType() + std::string(">");
  }

  static std::string EosPyType() { return std::string("MinimumEnergy") + T::EosPyType(); }

  MinimumEnergy() = default;

  auto GetOnDevice() { return MinimumEnergy<T>(t_.GetOnDevice()); }
  inline void Finalize() { t_.Finalize(); }

  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    const Real min_sie = t_.MinInternalEnergyFromDensity(rho);
    return t_.TemperatureFromDensityInternalEnergy(rho, std::max(sie, min_sie), lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temperature, Indexer_t &&lambda = nullptr) const {
    return t_.InternalEnergyFromDensityTemperature(rho, temperature, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real PressureFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    const Real min_sie = t_.MinInternalEnergyFromDensity(rho);
    return t_.PressureFromDensityInternalEnergy(rho, std::max(sie, min_sie), lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real
  MinInternalEnergyFromDensity(const Real rho, Indexer_t &&lambda = nullptr) const {
    return t_.MinInternalEnergyFromDensity(rho, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real EntropyFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    const Real min_sie = t_.MinInternalEnergyFromDensity(rho);
    return t_.EntropyFromDensityInternalEnergy(rho, std::max(sie, min_sie), lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    const Real min_sie = t_.MinInternalEnergyFromDensity(rho);
    return t_.SpecificHeatFromDensityInternalEnergy(rho, std::max(sie, min_sie), lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    const Real min_sie = t_.MinInternalEnergyFromDensity(rho);
    return t_.BulkModulusFromDensityInternalEnergy(rho, std::max(sie, min_sie), lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    const Real min_sie = t_.MinInternalEnergyFromDensity(rho);
    return t_.GruneisenParamFromDensityInternalEnergy(rho, std::max(sie, min_sie), lambda);
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
    if (output & thermalqs::specific_internal_energy) {
      t_.FillEos(rho, temp, energy, press, cv, bmod, output, lambda);
    } else {
      const Real min_sie = t_.MinInternalEnergyFromDensity(rho);
      t_.FillEos(rho, temp, std::max(energy, min_sie), press, cv, bmod, output, lambda);
    }
  }

  // vector implementations
  inline void choose_max_sie(const Real *sies, const Real *sie_use, const int num) const {
    // This code makes the assumption that `sie_use` is first populated with the
    // apprioriate minimum values
    static auto const name =
        singularity::mfuncname::member_func_name(typeid(MinimumEnergy<T>).name(), __func__);
    static auto const cname = name.c_str();
    portableFor(
        cname, 0, num,
        PORTABLE_LAMBDA(const int i) { sie_use[i] = std::max(sies[i], sie_use[i]); });
  }

  template <typename LambdaIndexer>
  inline void TemperatureFromDensityInternalEnergy(
      const Real *rhos, const Real *sies, Real *temperatures, Real *scratch,
      const int num, LambdaIndexer &&lambdas, Transform &&transform = Transform()) const {
    Real *sie_used = scratch;
    // First populate sies with minimum energies
    t_.MinInternalEnergyFromDensity(rhos, sie_used, num, lambdas);
    // Chose the maximum between the input and the minimum energy
    choose_max_sie(sies, sie_used, num);
    t_.TemperatureFromDensityInternalEnergy(
        rhos, sie_used, temperatures, &scratch[num], num,
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
    Real *sie_used = scratch;
    // First populate sies with minimum energies
    t_.MinInternalEnergyFromDensity(rhos, sie_used, num, lambdas);
    // Chose the maximum between the input and the minimum energy
    choose_max_sie(sies, sie_used, num);
    t_.PressureFromDensityInternalEnergy(rhos, sie_used, pressures, &scratch[num],
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
    Real *sie_used = scratch;
    // First populate sies with minimum energies
    t_.MinInternalEnergyFromDensity(rhos, sie_used, num, lambdas);
    // Chose the maximum between the input and the minimum energy
    choose_max_sie(sies, sie_used, num);
    t_.SpecificHeatFromDensityInternalEnergy(rhos, sie_used, cvs, &scratch[num], num,
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
    Real *sie_used = scratch;
    // First populate sies with minimum energies
    t_.MinInternalEnergyFromDensity(rhos, sie_used, num, lambdas);
    // Chose the maximum between the input and the minimum energy
    choose_max_sie(sies, sie_used, num);
    t_.BulkModulusFromDensityInternalEnergy(rhos, sie_used, bmods, &scratch[num], num,
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
    Real *sie_used = scratch;
    // First populate sies with minimum energies
    t_.MinInternalEnergyFromDensity(rhos, sie_used, num, lambdas);
    // Chose the maximum between the input and the minimum energy
    choose_max_sie(sies, sie_used, num);
    t_.GruneisenParamFromDensityInternalEnergy(rhos, sie_used, gm1s, &scratch[num],
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
    Real *sie_used = scratch;
    // First populate sies with minimum energies
    t_.MinInternalEnergyFromDensity(rhos, sie_used, num, lambdas);
    // Chose the maximum between the input and the minimum energy
    choose_max_sie(sies, sie_used, num);
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
  }
  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION void
  DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                       Indexer_t &&lambda, Real &rho, Real &sie) const {
    t_.DensityEnergyFromPressureTemperature(press, temp, lambda, rho, sie);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION void ValuesAtReferenceState(Real &rho, Real &temp, Real &sie,
                                                Real &press, Real &cv, Real &bmod,
                                                Real &dpde, Real &dvdt,
                                                Indexer_t &&lambda = nullptr) const {
    t_.ValuesAtReferenceState(rho, temp, sie, press, cv, bmod, dpde, dvdt, lambda);
  }

  PORTABLE_FORCEINLINE_FUNCTION Real MinimumDensity() const {
    return t_.MinimumDensity();
  }
  PORTABLE_FORCEINLINE_FUNCTION Real MinimumTemperature() const {
    return t_.MinimumTemperature();
  }

  inline constexpr bool IsModified() const { return true; }

  inline constexpr T UnmodifyOnce() { return t_; }

  inline constexpr decltype(auto) GetUnmodifiedObject() {
    return t_.GetUnmodifiedObject();
  }

 private:
  T t_;
};

} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_SHIFTED_EOS_

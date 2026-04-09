//------------------------------------------------------------------------------
// © 2021-2026 Triad National Security, LLC. All rights reserved.  This
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

// This file was created in part with generative AI

#ifndef _SINGULARITY_EOS_EOS_SHIFTED_EOS_
#define _SINGULARITY_EOS_EOS_SHIFTED_EOS_

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <utility>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_errors.hpp>
#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/eos_error.hpp>
#include <singularity-eos/eos/eos_base.hpp>
#include <singularity-eos/eos/modifiers/modifier_vector_macros.hpp>

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
  PORTABLE_FUNCTION void
  InternalEnergyFromDensityPressure(const Real rho, const Real P, Real &sie,
                                    Indexer_t &&lambda = nullptr) const {
    sie -= shift_;
    t_.InternalEnergyFromDensityPressure(rho, P, sie, lambda);
    sie += shift_;
  }

  template <typename Lambda_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  PTDerivativesFromPreferred(const Real rho, const Real sie, const Real P,
                             const Real temp, Lambda_t &&lambda, Real &dedP_T,
                             Real &drdP_T, Real &dedT_P, Real &drdT_P) const {
    t_.PTDerivativesFromPreferred(rho, sie - shift_, P, temp, lambda, dedP_T, drdP_T,
                                  dedT_P, drdT_P);
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
  template <typename Space, typename EnableIfSpace =
                                std::enable_if_t<!variadic_utils::has_int_index_v<Space>>>
  inline void shift_sies(const Space &s, const Real *sies, Real *shifted,
                         const int num) const {
    static auto const name =
        singularity::mfuncname::member_func_name(typeid(ShiftedEOS<T>).name(), __func__);
    static auto const cname = name.c_str();
    const auto shift_val = shift_;
    portableFor(
        cname, s, 0, num,
        PORTABLE_LAMBDA(const int i) { shifted[i] = sies[i] - shift_val; });
  }

  template <typename Space, typename EnableIfSpace =
                                std::enable_if_t<!variadic_utils::has_int_index_v<Space>>>
  inline void unshift_sies(const Space &s, Real *sies, const int num) const {
    static auto const name =
        singularity::mfuncname::member_func_name(typeid(ShiftedEOS<T>).name(), __func__);
    static auto const cname = name.c_str();
    const auto shift_val = shift_;
    portableFor(
        cname, s, 0, num, PORTABLE_LAMBDA(const int i) { sies[i] += shift_val; });
  }

  template <
      typename Space, typename LambdaIndexer,
      typename EnableIfSpace = std::enable_if_t<!variadic_utils::has_int_index_v<Space>>>
  inline void
  TemperatureFromDensityInternalEnergy(const Space &s, const Real *rhos, const Real *sies,
                                       Real *temperatures, Real *scratch, const int num,
                                       LambdaIndexer &&lambdas,
                                       Transform &&transform = Transform()) const {
    Real *shifted_sies = scratch;
    shift_sies(s, sies, shifted_sies, num);
    t_.TemperatureFromDensityInternalEnergy(
        s, rhos, shifted_sies, temperatures, &scratch[num], num,
        std::forward<LambdaIndexer>(lambdas), std::forward<Transform>(transform));
  }

  template <
      typename Space, typename LambdaIndexer,
      typename EnableIfSpace = std::enable_if_t<!variadic_utils::has_int_index_v<Space>>>
  inline void PressureFromDensityTemperature(const Space &s, const Real *rhos,
                                             const Real *temperatures, Real *pressures,
                                             Real *scratch, const int num,
                                             LambdaIndexer &&lambdas,
                                             Transform &&transform = Transform()) const {
    t_.PressureFromDensityTemperature(s, rhos, temperatures, pressures, scratch, num,
                                      std::forward<LambdaIndexer>(lambdas),
                                      std::forward<Transform>(transform));
  }

  template <
      typename Space, typename LambdaIndexer,
      typename EnableIfSpace = std::enable_if_t<!variadic_utils::has_int_index_v<Space>>>
  inline void PressureFromDensityInternalEnergy(
      const Space &s, const Real *rhos, const Real *sies, Real *pressures, Real *scratch,
      const int num, LambdaIndexer &&lambdas, Transform &&transform = Transform()) const {
    Real *shifted_sies = scratch;
    shift_sies(s, sies, shifted_sies, num);
    t_.PressureFromDensityInternalEnergy(s, rhos, shifted_sies, pressures, &scratch[num],
                                         num, std::forward<LambdaIndexer>(lambdas),
                                         std::forward<Transform>(transform));
  }

  template <
      typename Space, typename LambdaIndexer,
      typename EnableIfSpace = std::enable_if_t<!variadic_utils::has_int_index_v<Space>>>
  inline void MinInternalEnergyFromDensity(const Space &s, const Real *rhos, Real *sies,
                                           Real *scratch, const int num,
                                           LambdaIndexer &&lambdas,
                                           Transform &&transform = Transform()) const {
    t_.MinInternalEnergyFromDensity(s, rhos, sies, &scratch[num], num,
                                    std::forward<LambdaIndexer>(lambdas),
                                    std::forward<Transform>(transform));
    unshift_sies(s, sies, num);
  }

  template <
      typename Space, typename LambdaIndexer,
      typename EnableIfSpace = std::enable_if_t<!variadic_utils::has_int_index_v<Space>>>
  inline void
  SpecificHeatFromDensityTemperature(const Space &s, const Real *rhos,
                                     const Real *temperatures, Real *cvs, Real *scratch,
                                     const int num, LambdaIndexer &&lambdas,
                                     Transform &&transform = Transform()) const {
    t_.SpecificHeatFromDensityTemperature(s, rhos, temperatures, cvs, scratch, num,
                                          std::forward<LambdaIndexer>(lambdas),
                                          std::forward<Transform>(transform));
  }

  template <
      typename Space, typename LambdaIndexer,
      typename EnableIfSpace = std::enable_if_t<!variadic_utils::has_int_index_v<Space>>>
  inline void SpecificHeatFromDensityInternalEnergy(
      const Space &s, const Real *rhos, const Real *sies, Real *cvs, Real *scratch,
      const int num, LambdaIndexer &&lambdas, Transform &&transform = Transform()) const {
    Real *shifted_sies = scratch;
    shift_sies(s, sies, shifted_sies, num);
    t_.SpecificHeatFromDensityInternalEnergy(s, rhos, shifted_sies, cvs, &scratch[num],
                                             num, std::forward<LambdaIndexer>(lambdas),
                                             std::forward<Transform>(transform));
  }

  template <
      typename Space, typename LambdaIndexer,
      typename EnableIfSpace = std::enable_if_t<!variadic_utils::has_int_index_v<Space>>>
  inline void
  BulkModulusFromDensityTemperature(const Space &s, const Real *rhos,
                                    const Real *temperatures, Real *bmods, Real *scratch,
                                    const int num, LambdaIndexer &&lambdas,
                                    Transform &&transform = Transform()) const {
    t_.BulkModulusFromDensityTemperature(s, rhos, temperatures, bmods, scratch, num,
                                         std::forward<LambdaIndexer>(lambdas),
                                         std::forward<Transform>(transform));
  }

  template <
      typename Space, typename LambdaIndexer,
      typename EnableIfSpace = std::enable_if_t<!variadic_utils::has_int_index_v<Space>>>
  inline void BulkModulusFromDensityInternalEnergy(
      const Space &s, const Real *rhos, const Real *sies, Real *bmods, Real *scratch,
      const int num, LambdaIndexer &&lambdas, Transform &&transform = Transform()) const {
    Real *shifted_sies = scratch;
    shift_sies(s, sies, shifted_sies, num);
    t_.BulkModulusFromDensityInternalEnergy(s, rhos, shifted_sies, bmods, &scratch[num],
                                            num, std::forward<LambdaIndexer>(lambdas),
                                            std::forward<Transform>(transform));
  }

  template <
      typename Space, typename LambdaIndexer,
      typename EnableIfSpace = std::enable_if_t<!variadic_utils::has_int_index_v<Space>>>
  inline void GruneisenParamFromDensityTemperature(
      const Space &s, const Real *rhos, const Real *temperatures, Real *gm1s,
      Real *scratch, const int num, LambdaIndexer &&lambdas,
      Transform &&transform = Transform()) const {
    t_.GruneisenParamFromDensityTemperature(s, rhos, temperatures, gm1s, scratch, num,
                                            std::forward<LambdaIndexer>(lambdas),
                                            std::forward<Transform>(transform));
  }

  template <
      typename Space, typename LambdaIndexer,
      typename EnableIfSpace = std::enable_if_t<!variadic_utils::has_int_index_v<Space>>>
  inline void GruneisenParamFromDensityInternalEnergy(
      const Space &s, const Real *rhos, const Real *sies, Real *gm1s, Real *scratch,
      const int num, LambdaIndexer &&lambdas, Transform &&transform = Transform()) const {
    Real *shifted_sies = scratch;
    shift_sies(s, sies, shifted_sies, num);
    t_.GruneisenParamFromDensityInternalEnergy(s, rhos, shifted_sies, gm1s, &scratch[num],
                                               num, std::forward<LambdaIndexer>(lambdas),
                                               std::forward<Transform>(transform));
  }

  template <
      typename Space, typename LambdaIndexer,
      typename EnableIfSpace = std::enable_if_t<!variadic_utils::has_int_index_v<Space>>>
  inline void InternalEnergyFromDensityTemperature(
      const Space &s, const Real *rhos, const Real *temperatures, Real *sies,
      Real *scratch, const int num, LambdaIndexer &&lambdas,
      Transform &&transform = Transform()) const {
    t_.InternalEnergyFromDensityTemperature(s, rhos, temperatures, sies, scratch, num,
                                            std::forward<LambdaIndexer>(lambdas),
                                            std::forward<Transform>(transform));
    unshift_sies(s, sies, num);
  }

  template <
      typename Space, typename LambdaIndexer,
      typename EnableIfSpace = std::enable_if_t<!variadic_utils::has_int_index_v<Space>>>
  inline void EntropyFromDensityTemperature(const Space &s, const Real *rhos,
                                            const Real *temperatures, Real *entropies,
                                            Real *scratch, const int num,
                                            LambdaIndexer &&lambdas,
                                            Transform &&transform = Transform()) const {
    t_.EntropyFromDensityTemperature(s, rhos, temperatures, entropies, scratch, num,
                                     std::forward<LambdaIndexer>(lambdas),
                                     std::forward<Transform>(transform));
  }

  template <
      typename Space, typename LambdaIndexer,
      typename EnableIfSpace = std::enable_if_t<!variadic_utils::has_int_index_v<Space>>>
  inline void EntropyFromDensityInternalEnergy(
      const Space &s, const Real *rhos, const Real *sies, Real *entropies, Real *scratch,
      const int num, LambdaIndexer &&lambdas, Transform &&transform = Transform()) const {
    Real *shifted_sies = scratch;
    shift_sies(s, sies, shifted_sies, num);
    t_.EntropyFromDensityInternalEnergy(s, rhos, shifted_sies, entropies, &scratch[num],
                                        num, std::forward<LambdaIndexer>(lambdas),
                                        std::forward<Transform>(transform));
  }

  SG_MODIFIER_DEVICE_WRAP_ALL()

  constexpr static inline int nlambda() noexcept { return T::nlambda(); }
  template <typename Indexable>
  static inline constexpr bool NeedsLambda() {
    return T::template NeedsLambda<Indexable>();
  }

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

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real MeanAtomicMassFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return t_.MeanAtomicMassFromDensityTemperature(rho, temperature, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real MeanAtomicNumberFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return t_.MeanAtomicNumberFromDensityTemperature(rho, temperature, lambda);
  }

  SG_ADD_MODIFIER_METHODS(T, t_);
  SG_ADD_MODIFIER_MEAN_METHODS(t_);
  SG_ADD_MODIFIER_INTROSPECTION_METHODS(t_);

 private:
  T t_;
  double shift_;
};

} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_SHIFTED_EOS_

//------------------------------------------------------------------------------
// © 2021-2026. Triad National Security, LLC. All rights reserved.  This
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

// This file was made in part with generative AI.

#ifndef EOS_VARIANT_HPP
#define EOS_VARIANT_HPP

#include <memory>
#include <set>
#include <string>
#include <vector>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_errors.hpp>
#include <ports-of-call/variant.hpp>
#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/variadic_utils.hpp>
#include <singularity-eos/eos/eos_base.hpp>

using Real = double;

namespace singularity {

template <typename... Ts>
using eos_variant = PortsOfCall::variant<Ts...>;

// Provide default functionality when lambda isn't passed to vector
// functions For an example of what this might concretize to (albeit
// with fewer arguments than the functions targeted by this macro)
// look for the vector implementation of MinInternalEnergyFromDensity
// below.
// TODO(JMM): I decided to keep the names in the arguments to macro
// even though it's not strictly necessary, as I think it's more
// legible and produces more useful output in, e.g., a debugger.
struct NullIndexer {
  PORTABLE_FORCEINLINE_FUNCTION
  Real *operator[](int i) { return nullptr; }
  PORTABLE_FORCEINLINE_FUNCTION
  Real *operator[](int i) const { return nullptr; }
};

// Helper macros to reduce boilerplate in the Variant vector forwarding
// wrappers below.
#define SG_VARIANT_VEC_2IN_1OUT(NAME, IN1, IN2, OUT)                                     \
  template <typename Space, typename RealIndexer, typename ConstRealIndexer,             \
            typename = std::enable_if_t<!variadic_utils::has_int_index_v<Space>>>        \
  inline void NAME(const Space &s, ConstRealIndexer &&IN1, ConstRealIndexer &&IN2,       \
                   RealIndexer &&OUT, const int num) const {                             \
    NullIndexer lambdas{};                                                               \
    return NAME(s, std::forward<ConstRealIndexer>(IN1),                                  \
                std::forward<ConstRealIndexer>(IN2), std::forward<RealIndexer>(OUT),     \
                num, lambdas);                                                           \
  }                                                                                      \
  template <typename Space, typename RealIndexer, typename ConstRealIndexer,             \
            typename LambdaIndexer,                                                      \
            typename = std::enable_if_t<!variadic_utils::has_int_index_v<Space>>>        \
  inline void NAME(const Space &s, ConstRealIndexer &&IN1, ConstRealIndexer &&IN2,       \
                   RealIndexer &&OUT, const int num, LambdaIndexer &&lambdas) const {    \
    return PortsOfCall::visit(                                                           \
        [&s, &IN1, &IN2, &OUT, &num, &lambdas](const auto &eos) {                        \
          return eos.NAME(s, std::forward<ConstRealIndexer>(IN1),                        \
                          std::forward<ConstRealIndexer>(IN2),                           \
                          std::forward<RealIndexer>(OUT), num,                           \
                          std::forward<LambdaIndexer>(lambdas));                         \
        },                                                                               \
        eos_);                                                                           \
  }                                                                                      \
  template <typename Space, typename RealIndexer, typename ConstRealIndexer,             \
            typename = std::enable_if_t<!variadic_utils::has_int_index_v<Space>>>        \
  inline void NAME(const Space &s, ConstRealIndexer &&IN1, ConstRealIndexer &&IN2,       \
                   RealIndexer &&OUT, Real *scratch, const int num) const {              \
    NullIndexer lambdas{};                                                               \
    return NAME(s, std::forward<ConstRealIndexer>(IN1),                                  \
                std::forward<ConstRealIndexer>(IN2), std::forward<RealIndexer>(OUT),     \
                scratch, num, lambdas);                                                  \
  }                                                                                      \
  template <typename Space, typename RealIndexer, typename ConstRealIndexer,             \
            typename LambdaIndexer,                                                      \
            typename = std::enable_if_t<!variadic_utils::has_int_index_v<Space>>>        \
  inline void NAME(const Space &s, ConstRealIndexer &&IN1, ConstRealIndexer &&IN2,       \
                   RealIndexer &&OUT, Real *scratch, const int num,                      \
                   LambdaIndexer &&lambdas) const {                                      \
    return PortsOfCall::visit(                                                           \
        [&s, &IN1, &IN2, &OUT, &scratch, &num, &lambdas](const auto &eos) {              \
          return eos.NAME(s, std::forward<ConstRealIndexer>(IN1),                        \
                          std::forward<ConstRealIndexer>(IN2),                           \
                          std::forward<RealIndexer>(OUT), scratch, num,                  \
                          std::forward<LambdaIndexer>(lambdas));                         \
        },                                                                               \
        eos_);                                                                           \
  }                                                                                      \
  template <typename RealIndexer, typename ConstRealIndexer,                             \
            typename =                                                                   \
                std::enable_if_t<variadic_utils::has_int_index_v<ConstRealIndexer>>>     \
  inline void NAME(ConstRealIndexer &&IN1, ConstRealIndexer &&IN2, RealIndexer &&OUT,    \
                   const int num) const {                                                \
    return NAME(PortsOfCall::Exec::Device(), std::forward<ConstRealIndexer>(IN1),        \
                std::forward<ConstRealIndexer>(IN2), std::forward<RealIndexer>(OUT),     \
                num);                                                                    \
  }                                                                                      \
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer,     \
            typename =                                                                   \
                std::enable_if_t<variadic_utils::has_int_index_v<ConstRealIndexer>>>     \
  inline void NAME(ConstRealIndexer &&IN1, ConstRealIndexer &&IN2, RealIndexer &&OUT,    \
                   const int num, LambdaIndexer &&lambdas) const {                       \
    return NAME(PortsOfCall::Exec::Device(), std::forward<ConstRealIndexer>(IN1),        \
                std::forward<ConstRealIndexer>(IN2), std::forward<RealIndexer>(OUT),     \
                num, std::forward<LambdaIndexer>(lambdas));                              \
  }                                                                                      \
  template <typename RealIndexer, typename ConstRealIndexer,                             \
            typename =                                                                   \
                std::enable_if_t<variadic_utils::has_int_index_v<ConstRealIndexer>>>     \
  inline void NAME(ConstRealIndexer &&IN1, ConstRealIndexer &&IN2, RealIndexer &&OUT,    \
                   Real *scratch, const int num) const {                                 \
    return NAME(PortsOfCall::Exec::Device(), std::forward<ConstRealIndexer>(IN1),        \
                std::forward<ConstRealIndexer>(IN2), std::forward<RealIndexer>(OUT),     \
                scratch, num);                                                           \
  }                                                                                      \
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer,     \
            typename =                                                                   \
                std::enable_if_t<variadic_utils::has_int_index_v<ConstRealIndexer>>>     \
  inline void NAME(ConstRealIndexer &&IN1, ConstRealIndexer &&IN2, RealIndexer &&OUT,    \
                   Real *scratch, const int num, LambdaIndexer &&lambdas) const {        \
    return NAME(PortsOfCall::Exec::Device(), std::forward<ConstRealIndexer>(IN1),        \
                std::forward<ConstRealIndexer>(IN2), std::forward<RealIndexer>(OUT),     \
                scratch, num, std::forward<LambdaIndexer>(lambdas));                     \
  }
// This one does both FromDensityTemperature and
// FromDensityInternalEnergy at once. Not always useful, but
// frequently is.
#define SG_VARIANT_VEC_FOR(OUTNAME)                                                      \
  SG_VARIANT_VEC_2IN_1OUT(OUTNAME##FromDensityTemperature, rhos, temps, OUTNAME##s)      \
  SG_VARIANT_VEC_2IN_1OUT(OUTNAME##FromDensityInternalEnergy, rhos, sies, OUTNAME##s)

template <typename... EOSs>
class Variant {
 private:
  eos_variant<EOSs...> eos_;

 public:
  template <typename EOSChoice,
            typename std::enable_if<
                !std::is_same<Variant, typename std::decay<EOSChoice>::type>::value,
                bool>::type = true>
  PORTABLE_FUNCTION Variant(EOSChoice &&choice)
      : eos_(std::move(std::forward<EOSChoice>(choice))) {}

  Variant() = default;

  template <typename EOSChoice,
            typename std::enable_if<
                !std::is_same<Variant, typename std::decay<EOSChoice>::type>::value,
                bool>::type = true>
  PORTABLE_FUNCTION Variant &operator=(EOSChoice &&eos) {
    eos_ = std::move(std::forward<EOSChoice>(eos));
    return *this;
  }

  template <typename EOSChoice,
            typename std::enable_if<
                !std::is_same<Variant, typename std::decay<EOSChoice>::type>::value,
                bool>::type = true>
  PORTABLE_INLINE_FUNCTION EOSChoice get() const {
    return PortsOfCall::get<EOSChoice>(eos_);
  }

  Variant GetOnDevice() {
    return PortsOfCall::visit(
        [](auto &eos) { return eos_variant<EOSs...>(eos.GetOnDevice()); }, eos_);
  }

  // Place member functions here
  PORTABLE_INLINE_FUNCTION
  void CheckParams() const {
    return PortsOfCall::visit([](auto &eos) { return eos.CheckParams(); }, eos_);
  }

  template <typename Functor_t>
  PORTABLE_INLINE_FUNCTION void EvaluateDevice(const Functor_t f) const {
    return PortsOfCall::visit([&f](const auto &eos) { return eos.EvaluateDevice(f); },
                              eos_);
  }

  template <typename Functor_t>
  void EvaluateHost(Functor_t &f) const {
    return PortsOfCall::visit([&f](const auto &eos) { return eos.EvaluateHost(f); },
                              eos_);
  }

  // EOS modifier object-oriented API
  template <template <class> typename Mod>
  constexpr bool ModifiedInVariant() const {
    return PortsOfCall::visit(
        [](const auto &eos) { return eos.template ModifiedInList<Mod, EOSs...>(); },
        eos_);
  }
  template <template <class> typename Mod, typename... Args>
  constexpr auto Modify(Args &&...args) const {
    PORTABLE_ALWAYS_REQUIRE(ModifiedInVariant<Mod>(), "Modifier must be in variant");
    return PortsOfCall::visit(
        [&](const auto &eos) {
          auto modified = eos.template ConditionallyModify<Mod>(
              AvailableEOSTypes(), std::forward<Args>(args)...);
          return eos_variant<EOSs...>(modified);
        },
        eos_);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return PortsOfCall::visit(
        [&rho, &sie, &lambda](const auto &eos) {
          return eos.TemperatureFromDensityInternalEnergy(rho, sie, lambda);
        },
        eos_);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return PortsOfCall::visit(
        [&rho, &temperature, &lambda](const auto &eos) {
          return eos.InternalEnergyFromDensityTemperature(rho, temperature, lambda);
        },
        eos_);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return PortsOfCall::visit(
        [&rho, &temperature, &lambda](const auto &eos) {
          return eos.PressureFromDensityTemperature(rho, temperature, lambda);
        },
        eos_);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return PortsOfCall::visit(
        [&rho, &sie, &lambda](const auto &eos) {
          return eos.PressureFromDensityInternalEnergy(rho, sie, lambda);
        },
        eos_);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real MinInternalEnergyFromDensity(
      const Real rho, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return PortsOfCall::visit(
        [&rho, &lambda](const auto &eos) {
          return eos.MinInternalEnergyFromDensity(rho, lambda);
        },
        eos_);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real
  EntropyFromDensityTemperature(const Real rho, const Real temperature,
                                Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return PortsOfCall::visit(
        [&rho, &temperature, &lambda](const auto &eos) {
          return eos.EntropyFromDensityTemperature(rho, temperature, lambda);
        },
        eos_);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real EntropyFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return PortsOfCall::visit(
        [&rho, &sie, &lambda](const auto &eos) {
          return eos.EntropyFromDensityInternalEnergy(rho, sie, lambda);
        },
        eos_);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GibbsFreeEnergyFromDensityTemperature(
      const Real rho, const Real T,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return PortsOfCall::visit(
        [&rho, &T, &lambda](const auto &eos) {
          return eos.GibbsFreeEnergyFromDensityTemperature(rho, T, lambda);
        },
        eos_);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GibbsFreeEnergyFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return PortsOfCall::visit(
        [&rho, &sie, &lambda](const auto &eos) {
          return eos.GibbsFreeEnergyFromDensityInternalEnergy(rho, sie, lambda);
        },
        eos_);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return PortsOfCall::visit(
        [&rho, &temperature, &lambda](const auto &eos) {
          return eos.SpecificHeatFromDensityTemperature(rho, temperature, lambda);
        },
        eos_);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return PortsOfCall::visit(
        [&rho, &sie, &lambda](const auto &eos) {
          return eos.SpecificHeatFromDensityInternalEnergy(rho, sie, lambda);
        },
        eos_);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return PortsOfCall::visit(
        [&rho, &temperature, &lambda](const auto &eos) {
          return eos.BulkModulusFromDensityTemperature(rho, temperature, lambda);
        },
        eos_);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return PortsOfCall::visit(
        [&rho, &sie, &lambda](const auto &eos) {
          return eos.BulkModulusFromDensityInternalEnergy(rho, sie, lambda);
        },
        eos_);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return PortsOfCall::visit(
        [&rho, &temperature, &lambda](const auto &eos) {
          return eos.GruneisenParamFromDensityTemperature(rho, temperature, lambda);
        },
        eos_);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return PortsOfCall::visit(
        [&rho, &sie, &lambda](const auto &eos) {
          return eos.GruneisenParamFromDensityInternalEnergy(rho, sie, lambda);
        },
        eos_);
  }

  // TODO(JMM): Do we need a vectorized version of this call?
  template <typename Lambda_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  PTDerivativesFromPreferred(const Real rho, const Real sie, const Real P, const Real T,
                             Lambda_t &&lambda, Real &dedP_T, Real &drdP_T, Real &dedT_P,
                             Real &drdT_P) const {
    return PortsOfCall::visit(
        [&](const auto &eos) {
          return eos.PTDerivativesFromPreferred(rho, sie, P, T, lambda, dedP_T, drdP_T,
                                                dedT_P, drdT_P);
        },
        eos_);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  FillEos(Real &rho, Real &temp, Real &energy, Real &press, Real &cv, Real &bmod,
          const unsigned long output,
          Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return PortsOfCall::visit(
        [&rho, &temp, &energy, &press, &cv, &bmod, &output, &lambda](const auto &eos) {
          return eos.FillEos(rho, temp, energy, press, cv, bmod, output, lambda);
        },
        eos_);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  ReferenceDensityTemperature(Real &rho, Real &T,
                              Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return PortsOfCall::visit(
        [&rho, &T, &lambda](const auto &eos) {
          return eos.ReferenceDensityTemperature(rho, T, lambda);
        },
        eos_);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                         Real &bmod, Real &dpde, Real &dvdt,
                         Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return PortsOfCall::visit(
        [&rho, &temp, &sie, &press, &cv, &bmod, &dpde, &dvdt, &lambda](const auto &eos) {
          return eos.ValuesAtReferenceState(rho, temp, sie, press, cv, bmod, dpde, dvdt,
                                            lambda);
        },
        eos_);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                       Indexer_t &&lambda, Real &rho, Real &sie) const {
    return PortsOfCall::visit(
        [&press, &temp, &lambda, &rho, &sie](const auto &eos) {
          return eos.DensityEnergyFromPressureTemperature(press, temp, lambda, rho, sie);
        },
        eos_);
  }
  PORTABLE_INLINE_FUNCTION void DensityEnergyFromPressureTemperature(const Real press,
                                                                     const Real temp,
                                                                     Real &rho,
                                                                     Real &sie) const {
    return PortsOfCall::visit(
        [&press, &temp, &rho, &sie](const auto &eos) {
          return eos.DensityEnergyFromPressureTemperature(press, temp, rho, sie);
        },
        eos_);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void InternalEnergyFromDensityPressure(
      const Real rho, const Real P, Real &sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return PortsOfCall::visit(
        [&rho, &P, &sie, &lambda](const auto &eos) {
          return eos.InternalEnergyFromDensityPressure(rho, P, sie, lambda);
        },
        eos_);
  }

  PORTABLE_INLINE_FUNCTION
  Real RhoPmin(const Real temp) const {
    return PortsOfCall::visit([&temp](const auto &eos) { return eos.RhoPmin(temp); },
                              eos_);
  }

  PORTABLE_FORCEINLINE_FUNCTION
  Real MinimumDensity() const {
    return PortsOfCall::visit([](const auto &eos) { return eos.MinimumDensity(); }, eos_);
  }

  PORTABLE_FORCEINLINE_FUNCTION
  Real MinimumTemperature() const {
    return PortsOfCall::visit([](const auto &eos) { return eos.MinimumTemperature(); },
                              eos_);
  }

  PORTABLE_FORCEINLINE_FUNCTION
  Real MaximumDensity() const {
    return PortsOfCall::visit([](const auto &eos) { return eos.MaximumDensity(); }, eos_);
  }

  PORTABLE_FORCEINLINE_FUNCTION
  Real MinimumPressure() const {
    return PortsOfCall::visit([](const auto &eos) { return eos.MinimumPressure(); },
                              eos_);
  }

  PORTABLE_FORCEINLINE_FUNCTION
  Real MaximumPressureAtTemperature(const Real temp) const {
    return PortsOfCall::visit(
        [&temp](const auto &eos) { return eos.MaximumPressureAtTemperature(temp); },
        eos_);
  }

  // Atomic mass/atomic number functions
  PORTABLE_INLINE_FUNCTION
  Real MeanAtomicMass() const {
    return PortsOfCall::visit([](const auto &eos) { return eos.MeanAtomicMass(); }, eos_);
  }
  PORTABLE_INLINE_FUNCTION
  Real MeanAtomicNumber() const {
    return PortsOfCall::visit([](const auto &eos) { return eos.MeanAtomicNumber(); },
                              eos_);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real MeanAtomicMassFromDensityTemperature(
      const Real rho, const Real T,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return PortsOfCall::visit(
        [&rho, &T, &lambda](const auto &eos) {
          return eos.MeanAtomicMassFromDensityTemperature(rho, T, lambda);
        },
        eos_);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real MeanAtomicNumberFromDensityTemperature(
      const Real rho, const Real T,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return PortsOfCall::visit(
        [&rho, &T, &lambda](const auto &eos) {
          return eos.MeanAtomicNumberFromDensityTemperature(rho, T, lambda);
        },
        eos_);
  }

  /*
  Vector versions of the member functions run on the host but the scalar
  lookups will run on the device

  RealIndexer must have an operator[](int) that returns a Real. e.g., Real*
  ConstRealIndexer is as RealIndexer, but assumed const type.
  LambdaIndexer must have an operator[](int) that returns a Real*. e.g., Real**
  */
  SG_VARIANT_VEC_2IN_1OUT(TemperatureFromDensityInternalEnergy, rhos, sies, temperatures)
  SG_VARIANT_VEC_2IN_1OUT(InternalEnergyFromDensityTemperature, rhos, temperatures, sies)
  SG_VARIANT_VEC_FOR(Pressure)

  /// This is sort of what the SG_VARIANT_VEC would concretize too, though
  /// it has fewer arguments.
  template <typename Space, typename RealIndexer, typename ConstRealIndexer,
            typename = std::enable_if_t<!variadic_utils::has_int_index_v<Space>>>
  inline void MinInternalEnergyFromDensity(const Space &s, ConstRealIndexer &&rhos,
                                           RealIndexer &&sies, const int num) const {
    NullIndexer lambdas{}; // Returns null pointer for every index
    return MinInternalEnergyFromDensity(s, std::forward<ConstRealIndexer>(rhos),
                                        std::forward<RealIndexer>(sies), num, lambdas);
  }

  template <typename Space, typename RealIndexer, typename ConstRealIndexer,
            typename LambdaIndexer,
            typename = std::enable_if_t<!variadic_utils::has_int_index_v<Space>>>
  inline void MinInternalEnergyFromDensity(const Space &s, ConstRealIndexer &&rhos,
                                           RealIndexer &&sies, const int num,
                                           LambdaIndexer &&lambdas) const {
    return PortsOfCall::visit(
        [&s, &rhos, &sies, &num, &lambdas](const auto &eos) {
          return eos.MinInternalEnergyFromDensity(s, std::forward<ConstRealIndexer>(rhos),
                                                  std::forward<RealIndexer>(sies), num,
                                                  std::forward<LambdaIndexer>(lambdas));
        },
        eos_);
  }

  template <typename Space, typename RealIndexer, typename ConstRealIndexer,
            typename = std::enable_if_t<!variadic_utils::has_int_index_v<Space>>>
  inline void MinInternalEnergyFromDensity(const Space &s, ConstRealIndexer &&rhos,
                                           RealIndexer &&sies, Real *scratch,
                                           const int num) const {
    NullIndexer lambdas{}; // Returns null pointer for every index
    return MinInternalEnergyFromDensity(s, std::forward<ConstRealIndexer>(rhos),
                                        std::forward<RealIndexer>(sies), scratch, num,
                                        lambdas);
  }

  template <typename Space, typename RealIndexer, typename ConstRealIndexer,
            typename LambdaIndexer,
            typename = std::enable_if_t<!variadic_utils::has_int_index_v<Space>>>
  inline void MinInternalEnergyFromDensity(const Space &s, ConstRealIndexer &&rhos,
                                           RealIndexer &&sies, Real *scratch,
                                           const int num, LambdaIndexer &&lambdas) const {
    return PortsOfCall::visit(
        [&s, &rhos, &sies, &scratch, &num, &lambdas](const auto &eos) {
          return eos.MinInternalEnergyFromDensity(
              s, std::forward<ConstRealIndexer>(rhos), std::forward<RealIndexer>(sies),
              scratch, num, std::forward<LambdaIndexer>(lambdas));
        },
        eos_);
  }
  template <
      typename RealIndexer, typename ConstRealIndexer,
      typename = std::enable_if_t<variadic_utils::has_int_index_v<ConstRealIndexer>>>
  inline void MinInternalEnergyFromDensity(ConstRealIndexer &&rhos, RealIndexer &&sies,
                                           const int num) const {
    return MinInternalEnergyFromDensity(PortsOfCall::Exec::Device(),
                                        std::forward<ConstRealIndexer>(rhos),
                                        std::forward<RealIndexer>(sies), num);
  }

  template <
      typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer,
      typename = std::enable_if_t<variadic_utils::has_int_index_v<ConstRealIndexer>>>
  inline void MinInternalEnergyFromDensity(ConstRealIndexer &&rhos, RealIndexer &&sies,
                                           const int num, LambdaIndexer &&lambdas) const {
    return MinInternalEnergyFromDensity(
        PortsOfCall::Exec::Device(), std::forward<ConstRealIndexer>(rhos),
        std::forward<RealIndexer>(sies), num, std::forward<LambdaIndexer>(lambdas));
  }

  template <
      typename RealIndexer, typename ConstRealIndexer,
      typename = std::enable_if_t<variadic_utils::has_int_index_v<ConstRealIndexer>>>
  inline void MinInternalEnergyFromDensity(ConstRealIndexer &&rhos, RealIndexer &&sies,
                                           Real *scratch, const int num) const {
    return MinInternalEnergyFromDensity(PortsOfCall::Exec::Device(),
                                        std::forward<ConstRealIndexer>(rhos),
                                        std::forward<RealIndexer>(sies), scratch, num);
  }

  template <
      typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer,
      typename = std::enable_if_t<variadic_utils::has_int_index_v<ConstRealIndexer>>>
  inline void MinInternalEnergyFromDensity(ConstRealIndexer &&rhos, RealIndexer &&sies,
                                           Real *scratch, const int num,
                                           LambdaIndexer &&lambdas) const {
    return MinInternalEnergyFromDensity(PortsOfCall::Exec::Device(),
                                        std::forward<ConstRealIndexer>(rhos),
                                        std::forward<RealIndexer>(sies), scratch, num,
                                        std::forward<LambdaIndexer>(lambdas));
  }
  ///

  SG_VARIANT_VEC_FOR(Entropy)
  SG_VARIANT_VEC_FOR(GibbsFreeEnergy)
  SG_VARIANT_VEC_FOR(SpecificHeat)
  SG_VARIANT_VEC_FOR(BulkModulus)
  SG_VARIANT_VEC_FOR(GruneisenParam)
  SG_VARIANT_VEC_2IN_1OUT(InternalEnergyFromDensityPressure, rhos, Ps, sies)

  template <typename Space, typename RealIndexer,
            typename = std::enable_if_t<!variadic_utils::has_int_index_v<Space>>>
  inline void FillEos(const Space &s, RealIndexer &&rhos, RealIndexer &&temps,
                      RealIndexer &&energies, RealIndexer &&presses, RealIndexer &&cvs,
                      RealIndexer &&bmods, const int num,
                      const unsigned long output) const {
    NullIndexer lambdas{}; // Returns null pointer for every index
    return FillEos(s, std::forward<RealIndexer>(rhos), std::forward<RealIndexer>(temps),
                   std::forward<RealIndexer>(energies),
                   std::forward<RealIndexer>(presses), std::forward<RealIndexer>(cvs),
                   std::forward<RealIndexer>(bmods), num, output, lambdas);
  }

  template <typename Space, typename RealIndexer, typename LambdaIndexer,
            typename = std::enable_if_t<!variadic_utils::has_int_index_v<Space>>>
  inline void FillEos(const Space &s, RealIndexer &&rhos, RealIndexer &&temps,
                      RealIndexer &&energies, RealIndexer &&presses, RealIndexer &&cvs,
                      RealIndexer &&bmods, const int num, const unsigned long output,
                      LambdaIndexer &&lambdas) const {
    return PortsOfCall::visit(
        [&s, &rhos, &temps, &energies, &presses, &cvs, &bmods, &num, &output,
         &lambdas](const auto &eos) {
          return eos.FillEos(
              s, std::forward<RealIndexer>(rhos), std::forward<RealIndexer>(temps),
              std::forward<RealIndexer>(energies), std::forward<RealIndexer>(presses),
              std::forward<RealIndexer>(cvs), std::forward<RealIndexer>(bmods), num,
              output, std::forward<LambdaIndexer>(lambdas));
        },
        eos_);
  }

  template <typename RealIndexer,
            typename = std::enable_if_t<variadic_utils::has_int_index_v<RealIndexer>>>
  inline void FillEos(RealIndexer &&rhos, RealIndexer &&temps, RealIndexer &&energies,
                      RealIndexer &&presses, RealIndexer &&cvs, RealIndexer &&bmods,
                      const int num, const unsigned long output) const {
    return FillEos(PortsOfCall::Exec::Device(), std::forward<RealIndexer>(rhos),
                   std::forward<RealIndexer>(temps), std::forward<RealIndexer>(energies),
                   std::forward<RealIndexer>(presses), std::forward<RealIndexer>(cvs),
                   std::forward<RealIndexer>(bmods), num, output);
  }

  template <typename RealIndexer, typename LambdaIndexer,
            typename = std::enable_if_t<variadic_utils::has_int_index_v<RealIndexer>>>
  inline void FillEos(RealIndexer &&rhos, RealIndexer &&temps, RealIndexer &&energies,
                      RealIndexer &&presses, RealIndexer &&cvs, RealIndexer &&bmods,
                      const int num, const unsigned long output,
                      LambdaIndexer &&lambdas) const {
    return FillEos(PortsOfCall::Exec::Device(), std::forward<RealIndexer>(rhos),
                   std::forward<RealIndexer>(temps), std::forward<RealIndexer>(energies),
                   std::forward<RealIndexer>(presses), std::forward<RealIndexer>(cvs),
                   std::forward<RealIndexer>(bmods), num, output,
                   std::forward<LambdaIndexer>(lambdas));
  }

  // Serialization
  /*
    The methodology here is there are *three* size methods all EOS's provide:
    - `SharedMemorySizeInBytes()` which is the amount of memory a class can share
    - `DynamicMemorySizeInBytes()` which is the amount of memory not covered by
    `sizeof(this)`
    - `SerializedSizeInBytes()` which is the total size of the object.

    I wanted serialization machinery to work if you use a standalone
    class or if you use the variant. To make that possible, each class
    provides its own implementation of `SharedMemorySizeInBytes` and
    `DynamicMemorySizeInBytes()`. But then there is a separate
    implementation for the variant and for the base class for
    `SerializedSizeInBytes`, `Serialize`, and `DeSerialize`.
   */
  // JMM: This must be implemented separately for Variant vs the base
  // class/individual EOS's so that the variant state is properly
  // carried. Otherwise de-serialization would need to specify a type.
  std::size_t DynamicMemorySizeInBytes() const {
    return PortsOfCall::visit(
        [](const auto &eos) { return eos.DynamicMemorySizeInBytes(); }, eos_);
  }
  std::size_t DumpDynamicMemory(char *dst) {
    return PortsOfCall::visit([dst](auto &eos) { return eos.DumpDynamicMemory(dst); },
                              eos_);
  }
  std::size_t SetDynamicMemory(char *src,
                               const SharedMemSettings &stngs = DEFAULT_SHMEM_STNGS) {
    return PortsOfCall::visit(
        [src, stngs](auto &eos) { return eos.SetDynamicMemory(src, stngs); }, eos_);
  }
  std::size_t SharedMemorySizeInBytes() const {
    return PortsOfCall::visit(
        [](const auto &eos) { return eos.SharedMemorySizeInBytes(); }, eos_);
  }
  constexpr bool AllDynamicMemoryIsShareable() const {
    return PortsOfCall::visit(
        [](const auto &eos) { return eos.AllDynamicMemoryIsShareable(); }, eos_);
  }
  std::size_t SerializedSizeInBytes() const {
    return sizeof(*this) + DynamicMemorySizeInBytes();
  }
  std::size_t Serialize(char *dst) {
    memcpy(dst, this, sizeof(*this));
    std::size_t offst = sizeof(*this);
    std::size_t dyn_size = DynamicMemorySizeInBytes();
    if (dyn_size > 0) {
      offst += DumpDynamicMemory(dst + offst);
    }
    PORTABLE_ALWAYS_REQUIRE(offst == SerializedSizeInBytes(), "Serialization failed!");
    return offst;
  }
  auto Serialize() {
    std::size_t size = SerializedSizeInBytes();
    std::shared_ptr<char[]> dst(new char[size]);
    std::size_t new_size = Serialize(dst.get());
    PORTABLE_ALWAYS_REQUIRE(size == new_size, "Serialization failed!");
    return std::make_pair(size, dst);
  }
  std::size_t DeSerialize(char *src,
                          const SharedMemSettings &stngs = DEFAULT_SHMEM_STNGS) {
    memcpy(this, src, sizeof(*this));
    std::size_t offst = sizeof(*this);
    std::size_t dyn_size = DynamicMemorySizeInBytes();
    if (dyn_size > 0) {
      const bool sizes_same = AllDynamicMemoryIsShareable();
      if (stngs.CopyNeeded() && sizes_same) {
        memcpy(stngs.data, src + offst, dyn_size);
      }
      offst += SetDynamicMemory(src + offst, stngs);
    }
    return offst;
  }

  // Tooling for modifiers
  inline constexpr bool IsModified() const {
    return PortsOfCall::visit([](const auto &eos) { return eos.IsModified(); }, eos_);
  }

  inline constexpr Variant UnmodifyOnce() {
    return PortsOfCall::visit(
        [](auto &eos) -> eos_variant<EOSs...> {
          return eos_variant<EOSs...>(eos.UnmodifyOnce());
        },
        eos_);
  }

  inline constexpr Variant GetUnmodifiedObject() {
    return PortsOfCall::visit(
        [](auto &eos) -> eos_variant<EOSs...> {
          return eos_variant<EOSs...>(eos.GetUnmodifiedObject());
        },
        eos_);
  }

  PORTABLE_INLINE_FUNCTION
  unsigned long PreferredInput() const noexcept {
    return PortsOfCall::visit([](const auto &eos) { return eos.PreferredInput(); }, eos_);
  }

  PORTABLE_INLINE_FUNCTION
  unsigned long scratch_size(const std::string method, const unsigned int nelements) {
    return PortsOfCall::visit(
        [&](const auto &eos) { return eos.scratch_size(method, nelements); }, eos_);
  }

  PORTABLE_INLINE_FUNCTION
  unsigned long max_scratch_size(const unsigned int nelements) {
    return PortsOfCall::visit(
        [&](const auto &eos) { return eos.max_scratch_size(nelements); }, eos_);
  }

  PORTABLE_FORCEINLINE_FUNCTION
  int nlambda() noexcept {
    return PortsOfCall::visit([](const auto &eos) { return eos.nlambda(); }, eos_);
  }

  template <typename T>
  PORTABLE_INLINE_FUNCTION bool NeedsLambda() const {
    return PortsOfCall::visit(
        [](const auto &eos) { return eos.template NeedsLambda<T>(); }, eos_);
  }

  template <typename T>
  PORTABLE_INLINE_FUNCTION bool NeedsLambda(const T &t) const {
    return PortsOfCall::visit(
        [](const auto &eos) { return eos.template NeedsLambda<T>(); }, eos_);
  }

  template <typename T>
  PORTABLE_INLINE_FUNCTION bool IsType() const noexcept {
    return PortsOfCall::holds_alternative<T>(eos_);
  }

  PORTABLE_INLINE_FUNCTION
  void PrintParams() const noexcept {
    return PortsOfCall::visit([](const auto &eos) { return eos.PrintParams(); }, eos_);
  }

  inline std::string EosType() const {
    return PortsOfCall::visit([](const auto &eos) { return eos.EosType(); }, eos_);
  }

  inline std::string EosPyType() const {
    return PortsOfCall::visit([](const auto &eos) { return eos.EosPyType(); }, eos_);
  }

  template <typename Container_t = std::set<std::string>>
  static inline Container_t AvailableEOSs() {
    return {EOSs::EosType()...};
  }
  static inline constexpr auto AvailableEOSTypes() {
    return variadic_utils::type_list<EOSs...>();
  }

  inline void Finalize() noexcept {
    return PortsOfCall::visit([](auto &eos) { return eos.Finalize(); }, eos_);
  }
};
} // namespace singularity

#undef SG_VARIANT_VEC_FOR
#undef SG_VARIANT_VEC_2IN_1OUT
#endif // EOS_VARIANT_HPP

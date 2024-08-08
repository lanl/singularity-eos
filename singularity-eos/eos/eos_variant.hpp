//------------------------------------------------------------------------------
// © 2021-2023. Triad National Security, LLC. All rights reserved.  This
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

#ifndef EOS_VARIANT_HPP
#define EOS_VARIANT_HPP

#include <mpark/variant.hpp>
#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_errors.hpp>
#include <singularity-eos/base/variadic_utils.hpp>
#include <singularity-eos/eos/eos_base.hpp>

using Real = double;

namespace singularity {

template <typename... Ts>
using eos_variant = mpark::variant<Ts...>;

// Provide default functionality when lambda isn't passed to vector functions
struct NullIndexer {
  PORTABLE_FORCEINLINE_FUNCTION
  Real *operator[](int i) { return nullptr; }
  PORTABLE_FORCEINLINE_FUNCTION
  Real *operator[](int i) const { return nullptr; }
};

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

  Variant() noexcept = default;

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
    return mpark::get<EOSChoice>(eos_);
  }

  Variant GetOnDevice() {
    return mpark::visit([](auto &eos) { return eos_variant<EOSs...>(eos.GetOnDevice()); },
                        eos_);
  }

  // Place member functions here
  template <typename Functor_t>
  constexpr void Evaluate(Functor_t &f) const {
    return mpark::visit([&f](const auto &eos) { return eos.Evaluate(f); }, eos_);
  }

  // EOS modifier object-oriented API
  template <template <class> typename Mod>
  constexpr bool ModifiedInVariant() const {
    return mpark::visit(
        [](const auto &eos) { return eos.template ModifiedInList<Mod, EOSs...>(); },
        eos_);
  }
  template <template <class> typename Mod, typename... Args>
  constexpr auto Modify(Args &&...args) const {
    PORTABLE_ALWAYS_REQUIRE(ModifiedInVariant<Mod>(), "Modifier must be in variant");
    return mpark::visit(
        [&](const auto &eos) {
          auto modified = eos.template ConditionallyModify<Mod>(
              variadic_utils::type_list<EOSs...>(), std::forward<Args>(args)...);
          return eos_variant<EOSs...>(modified);
        },
        eos_);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return mpark::visit(
        [&rho, &sie, &lambda](const auto &eos) {
          return eos.TemperatureFromDensityInternalEnergy(rho, sie, lambda);
        },
        eos_);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return mpark::visit(
        [&rho, &temperature, &lambda](const auto &eos) {
          return eos.InternalEnergyFromDensityTemperature(rho, temperature, lambda);
        },
        eos_);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return mpark::visit(
        [&rho, &temperature, &lambda](const auto &eos) {
          return eos.PressureFromDensityTemperature(rho, temperature, lambda);
        },
        eos_);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return mpark::visit(
        [&rho, &sie, &lambda](const auto &eos) {
          return eos.PressureFromDensityInternalEnergy(rho, sie, lambda);
        },
        eos_);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real MinInternalEnergyFromDensity(
      const Real rho, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return mpark::visit(
        [&rho, &lambda](const auto &eos) {
          return eos.MinInternalEnergyFromDensity(rho, lambda);
        },
        eos_);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real
  EntropyFromDensityTemperature(const Real rho, const Real temperature,
                                Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return mpark::visit(
        [&rho, &temperature, &lambda](const auto &eos) {
          return eos.EntropyFromDensityTemperature(rho, temperature, lambda);
        },
        eos_);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real EntropyFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return mpark::visit(
        [&rho, &sie, &lambda](const auto &eos) {
          return eos.EntropyFromDensityInternalEnergy(rho, sie, lambda);
        },
        eos_);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return mpark::visit(
        [&rho, &temperature, &lambda](const auto &eos) {
          return eos.SpecificHeatFromDensityTemperature(rho, temperature, lambda);
        },
        eos_);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return mpark::visit(
        [&rho, &sie, &lambda](const auto &eos) {
          return eos.SpecificHeatFromDensityInternalEnergy(rho, sie, lambda);
        },
        eos_);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return mpark::visit(
        [&rho, &temperature, &lambda](const auto &eos) {
          return eos.BulkModulusFromDensityTemperature(rho, temperature, lambda);
        },
        eos_);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return mpark::visit(
        [&rho, &sie, &lambda](const auto &eos) {
          return eos.BulkModulusFromDensityInternalEnergy(rho, sie, lambda);
        },
        eos_);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return mpark::visit(
        [&rho, &temperature, &lambda](const auto &eos) {
          return eos.GruneisenParamFromDensityTemperature(rho, temperature, lambda);
        },
        eos_);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return mpark::visit(
        [&rho, &sie, &lambda](const auto &eos) {
          return eos.GruneisenParamFromDensityInternalEnergy(rho, sie, lambda);
        },
        eos_);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  FillEos(Real &rho, Real &temp, Real &energy, Real &press, Real &cv, Real &bmod,
          const unsigned long output,
          Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return mpark::visit(
        [&rho, &temp, &energy, &press, &cv, &bmod, &output, &lambda](const auto &eos) {
          return eos.FillEos(rho, temp, energy, press, cv, bmod, output, lambda);
        },
        eos_);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  ReferenceDensityTemperature(Real &rho, Real &T,
                              Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return mpark::visit(
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
    return mpark::visit(
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
    return mpark::visit(
        [&press, &temp, &lambda, &rho, &sie](const auto &eos) {
          return eos.DensityEnergyFromPressureTemperature(press, temp, lambda, rho, sie);
        },
        eos_);
  }

  PORTABLE_INLINE_FUNCTION
  Real RhoPmin(const Real temp) const {
    return mpark::visit([&temp](const auto &eos) { return eos.RhoPmin(temp); }, eos_);
  }

  PORTABLE_FORCEINLINE_FUNCTION
  Real MinimumDensity() const {
    return mpark::visit([](const auto &eos) { return eos.MinimumDensity(); }, eos_);
  }

  PORTABLE_FORCEINLINE_FUNCTION
  Real MinimumTemperature() const {
    return mpark::visit([](const auto &eos) { return eos.MinimumTemperature(); }, eos_);
  }

  /*
  Vector versions of the member functions run on the host but the scalar
  lookups will run on the device

  RealIndexer must have an operator[](int) that returns a Real. e.g., Real*
  ConstRealIndexer is as RealIndexer, but assumed const type.
  LambdaIndexer must have an operator[](int) that returns a Real*. e.g., Real**
  */
  template <typename RealIndexer, typename ConstRealIndexer>
  inline void
  TemperatureFromDensityInternalEnergy(ConstRealIndexer &&rhos, ConstRealIndexer &&sies,
                                       RealIndexer &&temperatures, const int num) const {
    NullIndexer lambdas{}; // Returns null pointer for every index
    return TemperatureFromDensityInternalEnergy(
        std::forward<ConstRealIndexer>(rhos), std::forward<ConstRealIndexer>(sies),
        std::forward<RealIndexer>(temperatures), num, lambdas);
  }

  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void
  TemperatureFromDensityInternalEnergy(ConstRealIndexer &&rhos, ConstRealIndexer &&sies,
                                       RealIndexer &&temperatures, const int num,
                                       LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &sies, &temperatures, &num, &lambdas](const auto &eos) {
          return eos.TemperatureFromDensityInternalEnergy(
              std::forward<ConstRealIndexer>(rhos), std::forward<ConstRealIndexer>(sies),
              std::forward<RealIndexer>(temperatures), num,
              std::forward<LambdaIndexer>(lambdas));
        },
        eos_);
  }
  template <typename RealIndexer, typename ConstRealIndexer>
  inline void TemperatureFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                                   ConstRealIndexer &&sies,
                                                   RealIndexer &&temperatures,
                                                   Real *scratch, const int num) const {
    NullIndexer lambdas{}; // Returns null pointer for every index
    return TemperatureFromDensityInternalEnergy(
        std::forward<ConstRealIndexer>(rhos), std::forward<ConstRealIndexer>(sies),
        std::forward<RealIndexer>(temperatures), scratch, num, lambdas);
  }

  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void
  TemperatureFromDensityInternalEnergy(ConstRealIndexer &&rhos, ConstRealIndexer &&sies,
                                       RealIndexer &&temperatures, Real *scratch,
                                       const int num, LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &sies, &temperatures, &scratch, &num, &lambdas](const auto &eos) {
          return eos.TemperatureFromDensityInternalEnergy(
              std::forward<ConstRealIndexer>(rhos), std::forward<ConstRealIndexer>(sies),
              std::forward<RealIndexer>(temperatures), scratch, num,
              std::forward<LambdaIndexer>(lambdas));
        },
        eos_);
  }

  template <typename RealIndexer, typename ConstRealIndexer>
  inline void InternalEnergyFromDensityTemperature(ConstRealIndexer &&rhos,
                                                   ConstRealIndexer &&temperatures,
                                                   RealIndexer &&sies,
                                                   const int num) const {
    NullIndexer lambdas{}; // Returns null pointer for every index
    return InternalEnergyFromDensityTemperature(
        std::forward<ConstRealIndexer>(rhos),
        std::forward<ConstRealIndexer>(temperatures), std::forward<RealIndexer>(sies),
        num, lambdas);
  }

  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void InternalEnergyFromDensityTemperature(ConstRealIndexer &&rhos,
                                                   ConstRealIndexer &&temperatures,
                                                   RealIndexer &&sies, const int num,
                                                   LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &temperatures, &sies, &num, &lambdas](const auto &eos) {
          return eos.InternalEnergyFromDensityTemperature(
              std::forward<ConstRealIndexer>(rhos),
              std::forward<ConstRealIndexer>(temperatures),
              std::forward<RealIndexer>(sies), num, std::forward<LambdaIndexer>(lambdas));
        },
        eos_);
  }

  template <typename RealIndexer, typename ConstRealIndexer>
  inline void InternalEnergyFromDensityTemperature(ConstRealIndexer &&rhos,
                                                   ConstRealIndexer &&temperatures,
                                                   RealIndexer &&sies, Real *scratch,
                                                   const int num) const {
    NullIndexer lambdas{}; // Returns null pointer for every index
    return InternalEnergyFromDensityTemperature(
        std::forward<ConstRealIndexer>(rhos),
        std::forward<ConstRealIndexer>(temperatures), sies, scratch, num, lambdas);
  }

  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void InternalEnergyFromDensityTemperature(ConstRealIndexer &&rhos,
                                                   ConstRealIndexer &&temperatures,
                                                   RealIndexer &&sies, Real *scratch,
                                                   const int num,
                                                   LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &temperatures, &sies, &scratch, &num, &lambdas](const auto &eos) {
          return eos.InternalEnergyFromDensityTemperature(
              std::forward<ConstRealIndexer>(rhos),
              std::forward<ConstRealIndexer>(temperatures), sies, scratch, num,
              std::forward<LambdaIndexer>(lambdas));
        },
        eos_);
  }

  template <typename RealIndexer, typename ConstRealIndexer>
  inline void
  PressureFromDensityTemperature(ConstRealIndexer &&rhos, ConstRealIndexer &&temperatures,
                                 RealIndexer &&pressures, const int num) const {
    NullIndexer lambdas{}; // Returns null pointer for every index
    return PressureFromDensityTemperature(std::forward<ConstRealIndexer>(rhos),
                                          std::forward<ConstRealIndexer>(temperatures),
                                          std::forward<RealIndexer>(pressures), num,
                                          lambdas);
  }

  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void PressureFromDensityTemperature(ConstRealIndexer &&rhos,
                                             ConstRealIndexer &&temperatures,
                                             RealIndexer &&pressures, const int num,
                                             LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &temperatures, &pressures, &num, &lambdas](const auto &eos) {
          return eos.PressureFromDensityTemperature(
              std::forward<ConstRealIndexer>(rhos),
              std::forward<ConstRealIndexer>(temperatures),
              std::forward<RealIndexer>(pressures), num,
              std::forward<LambdaIndexer>(lambdas));
        },
        eos_);
  }

  template <typename RealIndexer, typename ConstRealIndexer>
  inline void PressureFromDensityTemperature(ConstRealIndexer &&rhos,
                                             ConstRealIndexer &&temperatures,
                                             RealIndexer &&pressures, Real *scratch,
                                             const int num) const {
    NullIndexer lambdas{}; // Returns null pointer for every index
    return PressureFromDensityTemperature(std::forward<ConstRealIndexer>(rhos),
                                          std::forward<ConstRealIndexer>(temperatures),
                                          std::forward<RealIndexer>(pressures), scratch,
                                          num, lambdas);
  }

  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void
  PressureFromDensityTemperature(ConstRealIndexer &&rhos, ConstRealIndexer &&temperatures,
                                 RealIndexer &&pressures, Real *scratch, const int num,
                                 LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &temperatures, &pressures, &scratch, &num, &lambdas](const auto &eos) {
          return eos.PressureFromDensityTemperature(
              std::forward<ConstRealIndexer>(rhos),
              std::forward<ConstRealIndexer>(temperatures),
              std::forward<RealIndexer>(pressures), scratch, num,
              std::forward<LambdaIndexer>(lambdas));
        },
        eos_);
  }

  template <typename RealIndexer, typename ConstRealIndexer>
  inline void
  PressureFromDensityInternalEnergy(ConstRealIndexer &&rhos, ConstRealIndexer &&sies,
                                    RealIndexer &&pressures, const int num) const {
    NullIndexer lambdas{}; // Returns null pointer for every index
    return PressureFromDensityInternalEnergy(
        std::forward<ConstRealIndexer>(rhos), std::forward<ConstRealIndexer>(sies),
        std::forward<RealIndexer>(pressures), num, lambdas);
  }

  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void PressureFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                                ConstRealIndexer &&sies,
                                                RealIndexer &&pressures, const int num,
                                                LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &sies, &pressures, &num, &lambdas](const auto &eos) {
          return eos.PressureFromDensityInternalEnergy(
              std::forward<ConstRealIndexer>(rhos), std::forward<ConstRealIndexer>(sies),
              std::forward<RealIndexer>(pressures), num,
              std::forward<LambdaIndexer>(lambdas));
        },
        eos_);
  }

  template <typename RealIndexer, typename ConstRealIndexer>
  inline void PressureFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                                ConstRealIndexer &&sies,
                                                RealIndexer &&pressures, Real *scratch,
                                                const int num) const {
    NullIndexer lambdas{}; // Returns null pointer for every index
    return PressureFromDensityInternalEnergy(
        std::forward<ConstRealIndexer>(rhos), std::forward<ConstRealIndexer>(sies),
        std::forward<RealIndexer>(pressures), scratch, num, lambdas);
  }

  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void
  PressureFromDensityInternalEnergy(ConstRealIndexer &&rhos, ConstRealIndexer &&sies,
                                    RealIndexer &&pressures, Real *scratch, const int num,
                                    LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &sies, &pressures, &scratch, &num, &lambdas](const auto &eos) {
          return eos.PressureFromDensityInternalEnergy(
              std::forward<ConstRealIndexer>(rhos), std::forward<ConstRealIndexer>(sies),
              std::forward<RealIndexer>(pressures), scratch, num,
              std::forward<LambdaIndexer>(lambdas));
        },
        eos_);
  }
  ///
  template <typename RealIndexer, typename ConstRealIndexer>
  inline void MinInternalEnergyFromDensity(ConstRealIndexer &&rhos, RealIndexer &&sies,
                                           const int num) const {
    NullIndexer lambdas{}; // Returns null pointer for every index
    return MinInternalEnergyFromDensity(std::forward<ConstRealIndexer>(rhos),
                                        std::forward<RealIndexer>(sies), num, lambdas);
  }

  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void MinInternalEnergyFromDensity(ConstRealIndexer &&rhos, RealIndexer &&sies,
                                           const int num, LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &sies, &num, &lambdas](const auto &eos) {
          return eos.MinInternalEnergyFromDensity(std::forward<ConstRealIndexer>(rhos),
                                                  std::forward<RealIndexer>(sies), num,
                                                  std::forward<LambdaIndexer>(lambdas));
        },
        eos_);
  }

  template <typename RealIndexer, typename ConstRealIndexer>
  inline void MinInternalEnergyFromDensity(ConstRealIndexer &&rhos, RealIndexer &&sies,
                                           Real *scratch, const int num) const {
    NullIndexer lambdas{}; // Returns null pointer for every index
    return MinInternalEnergyFromDensity(std::forward<ConstRealIndexer>(rhos),
                                        std::forward<RealIndexer>(sies), scratch, num,
                                        lambdas);
  }

  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void MinInternalEnergyFromDensity(ConstRealIndexer &&rhos, RealIndexer &&sies,
                                           Real *scratch, const int num,
                                           LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &sies, &scratch, &num, &lambdas](const auto &eos) {
          return eos.MinInternalEnergyFromDensity(
              std::forward<ConstRealIndexer>(rhos), std::forward<RealIndexer>(sies),
              scratch, num, std::forward<LambdaIndexer>(lambdas));
        },
        eos_);
  }
  ///
  template <typename RealIndexer, typename ConstRealIndexer>
  inline void
  EntropyFromDensityTemperature(ConstRealIndexer &&rhos, ConstRealIndexer &&temperatures,
                                RealIndexer &&entropies, const int num) const {
    NullIndexer lambdas{}; // Returns null pointer for every index
    return EntropyFromDensityTemperature(std::forward<ConstRealIndexer>(rhos),
                                         std::forward<ConstRealIndexer>(temperatures),
                                         std::forward<RealIndexer>(entropies), num,
                                         lambdas);
  }

  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void EntropyFromDensityTemperature(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&temperatures,
                                            RealIndexer &&entropies, const int num,
                                            LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &temperatures, &entropies, &num, &lambdas](const auto &eos) {
          return eos.EntropyFromDensityTemperature(
              std::forward<ConstRealIndexer>(rhos),
              std::forward<ConstRealIndexer>(temperatures),
              std::forward<RealIndexer>(entropies), num,
              std::forward<LambdaIndexer>(lambdas));
        },
        eos_);
  }

  template <typename RealIndexer, typename ConstRealIndexer>
  inline void EntropyFromDensityTemperature(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&temperatures,
                                            RealIndexer &&entropies, Real *scratch,
                                            const int num) const {
    NullIndexer lambdas{}; // Returns null pointer for every index
    return EntropyFromDensityTemperature(std::forward<ConstRealIndexer>(rhos),
                                         std::forward<ConstRealIndexer>(temperatures),
                                         std::forward<RealIndexer>(entropies), scratch,
                                         num, lambdas);
  }

  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void
  EntropyFromDensityTemperature(ConstRealIndexer &&rhos, ConstRealIndexer &&temperatures,
                                RealIndexer &&entropies, Real *scratch, const int num,
                                LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &temperatures, &entropies, &scratch, &num, &lambdas](const auto &eos) {
          return eos.EntropyFromDensityTemperature(
              std::forward<ConstRealIndexer>(rhos),
              std::forward<ConstRealIndexer>(temperatures),
              std::forward<RealIndexer>(entropies), scratch, num,
              std::forward<LambdaIndexer>(lambdas));
        },
        eos_);
  }

  template <typename RealIndexer, typename ConstRealIndexer>
  inline void
  EntropyFromDensityInternalEnergy(ConstRealIndexer &&rhos, ConstRealIndexer &&sies,
                                   RealIndexer &&entropies, const int num) const {
    NullIndexer lambdas{}; // Returns null pointer for every index
    return EntropyFromDensityInternalEnergy(
        std::forward<ConstRealIndexer>(rhos), std::forward<ConstRealIndexer>(sies),
        std::forward<RealIndexer>(entropies), num, lambdas);
  }

  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void EntropyFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                               ConstRealIndexer &&sies,
                                               RealIndexer &&entropies, const int num,
                                               LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &sies, &entropies, &num, &lambdas](const auto &eos) {
          return eos.EntropyFromDensityInternalEnergy(
              std::forward<ConstRealIndexer>(rhos), std::forward<ConstRealIndexer>(sies),
              std::forward<RealIndexer>(entropies), num,
              std::forward<LambdaIndexer>(lambdas));
        },
        eos_);
  }

  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void EntropyFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                               ConstRealIndexer &&sies,
                                               RealIndexer &&entropies, Real *scratch,
                                               const int num) const {
    NullIndexer lambdas{}; // Returns null pointer for every index
    return EntropyFromDensityInternalEnergy(
        std::forward<ConstRealIndexer>(rhos), std::forward<ConstRealIndexer>(sies),
        std::forward<RealIndexer>(entropies), scratch, num, lambdas);
  }

  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void
  EntropyFromDensityInternalEnergy(ConstRealIndexer &&rhos, ConstRealIndexer &&sies,
                                   RealIndexer &&entropies, Real *scratch, const int num,
                                   LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &sies, &entropies, &scratch, &num, &lambdas](const auto &eos) {
          return eos.EntropyFromDensityInternalEnergy(
              std::forward<ConstRealIndexer>(rhos), std::forward<ConstRealIndexer>(sies),
              std::forward<RealIndexer>(entropies), scratch, num,
              std::forward<LambdaIndexer>(lambdas));
        },
        eos_);
  }

  template <typename RealIndexer, typename ConstRealIndexer>
  inline void SpecificHeatFromDensityTemperature(ConstRealIndexer &&rhos,
                                                 ConstRealIndexer &&temperatures,
                                                 RealIndexer &&cvs, const int num) const {
    NullIndexer lambdas{}; // Returns null pointer for every index
    return SpecificHeatFromDensityTemperature(
        std::forward<ConstRealIndexer>(rhos),
        std::forward<ConstRealIndexer>(temperatures), std::forward<RealIndexer>(cvs), num,
        lambdas);
  }

  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void SpecificHeatFromDensityTemperature(ConstRealIndexer &&rhos,
                                                 ConstRealIndexer &&temperatures,
                                                 RealIndexer &&cvs, const int num,
                                                 LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &temperatures, &cvs, &num, &lambdas](const auto &eos) {
          return eos.SpecificHeatFromDensityTemperature(
              std::forward<ConstRealIndexer>(rhos),
              std::forward<ConstRealIndexer>(temperatures),
              std::forward<RealIndexer>(cvs), num, std::forward<LambdaIndexer>(lambdas));
        },
        eos_);
  }

  template <typename RealIndexer, typename ConstRealIndexer>
  inline void SpecificHeatFromDensityTemperature(ConstRealIndexer &&rhos,
                                                 ConstRealIndexer &&temperatures,
                                                 RealIndexer &&cvs, Real *scratch,
                                                 const int num) const {
    NullIndexer lambdas{}; // Returns null pointer for every index
    return SpecificHeatFromDensityTemperature(
        std::forward<ConstRealIndexer>(rhos),
        std::forward<ConstRealIndexer>(temperatures), std::forward<RealIndexer>(cvs),
        scratch, num, lambdas);
  }

  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void SpecificHeatFromDensityTemperature(ConstRealIndexer &&rhos,
                                                 ConstRealIndexer &&temperatures,
                                                 RealIndexer &&cvs, Real *scratch,
                                                 const int num,
                                                 LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &temperatures, &cvs, &scratch, &num, &lambdas](const auto &eos) {
          return eos.SpecificHeatFromDensityTemperature(
              std::forward<ConstRealIndexer>(rhos),
              std::forward<ConstRealIndexer>(temperatures),
              std::forward<RealIndexer>(cvs), scratch, num,
              std::forward<LambdaIndexer>(lambdas));
        },
        eos_);
  }

  template <typename RealIndexer, typename ConstRealIndexer>
  inline void
  SpecificHeatFromDensityInternalEnergy(ConstRealIndexer &&rhos, ConstRealIndexer &&sies,
                                        RealIndexer &&cvs, const int num) const {
    NullIndexer lambdas{}; // Returns null pointer for every index
    return SpecificHeatFromDensityInternalEnergy(
        std::forward<ConstRealIndexer>(rhos), std::forward<ConstRealIndexer>(sies),
        std::forward<RealIndexer>(cvs), num, lambdas);
  }

  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void SpecificHeatFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                                    ConstRealIndexer &&sies,
                                                    RealIndexer &&cvs, const int num,
                                                    LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &sies, &cvs, &num, &lambdas](const auto &eos) {
          return eos.SpecificHeatFromDensityInternalEnergy(
              std::forward<ConstRealIndexer>(rhos), std::forward<ConstRealIndexer>(sies),
              std::forward<RealIndexer>(cvs), num, std::forward<LambdaIndexer>(lambdas));
        },
        eos_);
  }
  template <typename RealIndexer, typename ConstRealIndexer>
  inline void SpecificHeatFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                                    ConstRealIndexer &&sies,
                                                    RealIndexer &&cvs, Real *scratch,
                                                    const int num) const {
    NullIndexer lambdas{}; // Returns null pointer for every index
    return SpecificHeatFromDensityInternalEnergy(
        std::forward<ConstRealIndexer>(rhos), std::forward<ConstRealIndexer>(sies),
        std::forward<RealIndexer>(cvs), scratch, num, lambdas);
  }

  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void
  SpecificHeatFromDensityInternalEnergy(ConstRealIndexer &&rhos, ConstRealIndexer &&sies,
                                        RealIndexer &&cvs, Real *scratch, const int num,
                                        LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &sies, &cvs, &scratch, &num, &lambdas](const auto &eos) {
          return eos.SpecificHeatFromDensityInternalEnergy(
              std::forward<ConstRealIndexer>(rhos), std::forward<ConstRealIndexer>(sies),
              std::forward<RealIndexer>(cvs), scratch, num,
              std::forward<LambdaIndexer>(lambdas));
        },
        eos_);
  }

  template <typename RealIndexer, typename ConstRealIndexer>
  inline void BulkModulusFromDensityTemperature(ConstRealIndexer &&rhos,
                                                ConstRealIndexer &&temperatures,
                                                RealIndexer &&bmods,
                                                const int num) const {
    NullIndexer lambdas{}; // Returns null pointer for every index
    return BulkModulusFromDensityTemperature(std::forward<ConstRealIndexer>(rhos),
                                             std::forward<ConstRealIndexer>(temperatures),
                                             std::forward<RealIndexer>(bmods), num,
                                             lambdas);
  }

  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void BulkModulusFromDensityTemperature(ConstRealIndexer &&rhos,
                                                ConstRealIndexer &&temperatures,
                                                RealIndexer &&bmods, const int num,
                                                LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &temperatures, &bmods, &num, &lambdas](const auto &eos) {
          return eos.BulkModulusFromDensityTemperature(
              std::forward<ConstRealIndexer>(rhos),
              std::forward<ConstRealIndexer>(temperatures),
              std::forward<RealIndexer>(bmods), num,
              std::forward<LambdaIndexer>(lambdas));
        },
        eos_);
  }
  template <typename RealIndexer, typename ConstRealIndexer>
  inline void BulkModulusFromDensityTemperature(ConstRealIndexer &&rhos,
                                                ConstRealIndexer &&temperatures,
                                                RealIndexer &&bmods, Real *scratch,
                                                const int num) const {
    NullIndexer lambdas{}; // Returns null pointer for every index
    return BulkModulusFromDensityTemperature(std::forward<ConstRealIndexer>(rhos),
                                             std::forward<ConstRealIndexer>(temperatures),
                                             std::forward<RealIndexer>(bmods), scratch,
                                             num, lambdas);
  }

  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void BulkModulusFromDensityTemperature(ConstRealIndexer &&rhos,
                                                ConstRealIndexer &&temperatures,
                                                RealIndexer &&bmods, Real *scratch,
                                                const int num,
                                                LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &temperatures, &bmods, &scratch, &num, &lambdas](const auto &eos) {
          return eos.BulkModulusFromDensityTemperature(
              std::forward<ConstRealIndexer>(rhos),
              std::forward<ConstRealIndexer>(temperatures),
              std::forward<RealIndexer>(bmods), scratch, num,
              std::forward<LambdaIndexer>(lambdas));
        },
        eos_);
  }

  template <typename RealIndexer, typename ConstRealIndexer>
  inline void
  BulkModulusFromDensityInternalEnergy(ConstRealIndexer &&rhos, ConstRealIndexer &&sies,
                                       RealIndexer &&bmods, const int num) const {
    NullIndexer lambdas{}; // Returns null pointer for every index
    return BulkModulusFromDensityInternalEnergy(
        std::forward<ConstRealIndexer>(rhos), std::forward<ConstRealIndexer>(sies),
        std::forward<RealIndexer>(bmods), num, lambdas);
  }

  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void BulkModulusFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                                   ConstRealIndexer &&sies,
                                                   RealIndexer &&bmods, const int num,
                                                   LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &sies, &bmods, &num, &lambdas](const auto &eos) {
          return eos.BulkModulusFromDensityInternalEnergy(
              std::forward<ConstRealIndexer>(rhos), std::forward<ConstRealIndexer>(sies),
              std::forward<RealIndexer>(bmods), num,
              std::forward<LambdaIndexer>(lambdas));
        },
        eos_);
  }

  template <typename RealIndexer, typename ConstRealIndexer>
  inline void BulkModulusFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                                   ConstRealIndexer &&sies,
                                                   RealIndexer &&bmods, Real *scratch,
                                                   const int num) const {
    NullIndexer lambdas{}; // Returns null pointer for every index
    return BulkModulusFromDensityInternalEnergy(
        std::forward<ConstRealIndexer>(rhos), std::forward<ConstRealIndexer>(sies),
        std::forward<RealIndexer>(bmods), scratch, num, lambdas);
  }

  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void
  BulkModulusFromDensityInternalEnergy(ConstRealIndexer &&rhos, ConstRealIndexer &&sies,
                                       RealIndexer &&bmods, Real *scratch, const int num,
                                       LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &sies, &bmods, &scratch, &num, &lambdas](const auto &eos) {
          return eos.BulkModulusFromDensityInternalEnergy(
              std::forward<ConstRealIndexer>(rhos), std::forward<ConstRealIndexer>(sies),
              std::forward<RealIndexer>(bmods), scratch, num,
              std::forward<LambdaIndexer>(lambdas));
        },
        eos_);
  }

  template <typename RealIndexer, typename ConstRealIndexer>
  inline void GruneisenParamFromDensityTemperature(ConstRealIndexer &&rhos,
                                                   ConstRealIndexer &&temperatures,
                                                   RealIndexer &&gm1s,
                                                   const int num) const {
    NullIndexer lambdas{}; // Returns null pointer for every index
    return GruneisenParamFromDensityTemperature(
        std::forward<ConstRealIndexer>(rhos),
        std::forward<ConstRealIndexer>(temperatures), std::forward<RealIndexer>(gm1s),
        num, lambdas);
  }

  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void GruneisenParamFromDensityTemperature(ConstRealIndexer &&rhos,
                                                   ConstRealIndexer &&temperatures,
                                                   RealIndexer &&gm1s, const int num,
                                                   LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &temperatures, &gm1s, &num, &lambdas](const auto &eos) {
          return eos.GruneisenParamFromDensityTemperature(
              std::forward<ConstRealIndexer>(rhos),
              std::forward<ConstRealIndexer>(temperatures),
              std::forward<RealIndexer>(gm1s), num, std::forward<LambdaIndexer>(lambdas));
        },
        eos_);
  }

  template <typename RealIndexer, typename ConstRealIndexer>
  inline void GruneisenParamFromDensityTemperature(ConstRealIndexer &&rhos,
                                                   ConstRealIndexer &&temperatures,
                                                   RealIndexer &&gm1s, Real *scratch,
                                                   const int num) const {
    NullIndexer lambdas{}; // Returns null pointer for every index
    return GruneisenParamFromDensityTemperature(
        std::forward<ConstRealIndexer>(rhos),
        std::forward<ConstRealIndexer>(temperatures), std::forward<RealIndexer>(gm1s),
        scratch, num, lambdas);
  }

  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void GruneisenParamFromDensityTemperature(ConstRealIndexer &&rhos,
                                                   ConstRealIndexer &&temperatures,
                                                   RealIndexer &&gm1s, Real *scratch,
                                                   const int num,
                                                   LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &temperatures, &gm1s, &scratch, &num, &lambdas](const auto &eos) {
          return eos.GruneisenParamFromDensityTemperature(
              std::forward<ConstRealIndexer>(rhos),
              std::forward<ConstRealIndexer>(temperatures),
              std::forward<RealIndexer>(gm1s), scratch, num,
              std::forward<LambdaIndexer>(lambdas));
        },
        eos_);
  }

  template <typename RealIndexer, typename ConstRealIndexer>
  inline void GruneisenParamFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                                      ConstRealIndexer &&sies,
                                                      RealIndexer &&gm1s,
                                                      const int num) const {
    NullIndexer lambdas{}; // Returns null pointer for every index
    return GruneisenParamFromDensityInternalEnergy(std::forward<ConstRealIndexer>(rhos),
                                                   std::forward<ConstRealIndexer>(sies),
                                                   gm1s, num, lambdas);
  }

  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void GruneisenParamFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                                      ConstRealIndexer &&sies,
                                                      RealIndexer &&gm1s, const int num,
                                                      LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &sies, &gm1s, &lambdas, &num](const auto &eos) {
          return eos.GruneisenParamFromDensityInternalEnergy(
              std::forward<ConstRealIndexer>(rhos), std::forward<ConstRealIndexer>(sies),
              std::forward<RealIndexer>(gm1s), num, std::forward<LambdaIndexer>(lambdas));
        },
        eos_);
  }

  template <typename RealIndexer, typename ConstRealIndexer>
  inline void GruneisenParamFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                                      ConstRealIndexer &&sies,
                                                      RealIndexer &&gm1s, Real *scratch,
                                                      const int num) const {
    NullIndexer lambdas{}; // Returns null pointer for every index
    return GruneisenParamFromDensityInternalEnergy(
        std::forward<ConstRealIndexer>(rhos), std::forward<ConstRealIndexer>(sies),
        std::forward<RealIndexer>(gm1s), scratch, num, lambdas);
  }

  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void GruneisenParamFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                                      ConstRealIndexer &&sies,
                                                      RealIndexer &&gm1s, Real *scratch,
                                                      const int num,
                                                      LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &sies, &gm1s, &scratch, &lambdas, &num](const auto &eos) {
          return eos.GruneisenParamFromDensityInternalEnergy(
              std::forward<ConstRealIndexer>(rhos), std::forward<ConstRealIndexer>(sies),
              std::forward<RealIndexer>(gm1s), scratch, num,
              std::forward<LambdaIndexer>(lambdas));
        },
        eos_);
  }

  template <typename RealIndexer>
  inline void FillEos(RealIndexer &&rhos, RealIndexer &&temps, RealIndexer &&energies,
                      RealIndexer &&presses, RealIndexer &&cvs, RealIndexer &&bmods,
                      const int num, const unsigned long output) const {
    NullIndexer lambdas{}; // Returns null pointer for every index
    return FillEos(std::forward<RealIndexer>(rhos), std::forward<RealIndexer>(temps),
                   std::forward<RealIndexer>(energies),
                   std::forward<RealIndexer>(presses), std::forward<RealIndexer>(cvs),
                   std::forward<RealIndexer>(bmods), num, output, lambdas);
  }

  template <typename RealIndexer, typename LambdaIndexer>
  inline void FillEos(RealIndexer &&rhos, RealIndexer &&temps, RealIndexer &&energies,
                      RealIndexer &&presses, RealIndexer &&cvs, RealIndexer &&bmods,
                      const int num, const unsigned long output,
                      LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &temps, &energies, &presses, &cvs, &bmods, &num, &output,
         &lambdas](const auto &eos) {
          return eos.FillEos(
              std::forward<RealIndexer>(rhos), std::forward<RealIndexer>(temps),
              std::forward<RealIndexer>(energies), std::forward<RealIndexer>(presses),
              std::forward<RealIndexer>(cvs), std::forward<RealIndexer>(bmods), num,
              output, std::forward<LambdaIndexer>(lambdas));
        },
        eos_);
  }

  // Tooling for modifiers
  inline constexpr bool IsModified() const {
    return mpark::visit([](const auto &eos) { return eos.IsModified(); }, eos_);
  }

  inline constexpr Variant UnmodifyOnce() {
    return mpark::visit(
        [](auto &eos) -> eos_variant<EOSs...> {
          return eos_variant<EOSs...>(eos.UnmodifyOnce());
        },
        eos_);
  }

  inline constexpr Variant GetUnmodifiedObject() {
    return mpark::visit(
        [](auto &eos) -> eos_variant<EOSs...> {
          return eos_variant<EOSs...>(eos.GetUnmodifiedObject());
        },
        eos_);
  }

  PORTABLE_INLINE_FUNCTION
  unsigned long PreferredInput() const noexcept {
    return mpark::visit([](const auto &eos) { return eos.PreferredInput(); }, eos_);
  }

  PORTABLE_INLINE_FUNCTION
  unsigned long scratch_size(std::string method, unsigned int nelements) {
    return mpark::visit(
        [&](const auto &eos) { return eos.scratch_size(method, nelements); }, eos_);
  }

  PORTABLE_INLINE_FUNCTION
  unsigned long max_scratch_size(unsigned int nelements) {
    return mpark::visit([&](const auto &eos) { return eos.max_scratch_size(nelements); },
                        eos_);
  }

  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept {
    return mpark::visit([](const auto &eos) { return eos.nlambda(); }, eos_);
  }

  template <typename T>
  PORTABLE_INLINE_FUNCTION bool IsType() const noexcept {
    return mpark::holds_alternative<T>(eos_);
  }

  PORTABLE_INLINE_FUNCTION
  void PrintParams() const noexcept {
    return mpark::visit([](const auto &eos) { return eos.PrintParams(); }, eos_);
  }

  inline void Finalize() noexcept {
    return mpark::visit([](auto &eos) { return eos.Finalize(); }, eos_);
  }
};
} // namespace singularity

#endif // EOS_VARIANT_HPP

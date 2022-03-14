//------------------------------------------------------------------------------
// © 2021. Triad National Security, LLC. All rights reserved.  This
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

#include <ports-of-call/portability.hpp>
#include <variant/include/mpark/variant.hpp>
#include <singularity-eos/eos/eos_base.hpp>

using Real = double;

namespace singularity {

template <typename... Ts>
using eos_variant = mpark::variant<Ts...>;

// Provide default functionality when lambda isn't passed to vector functions
class NullIndexer {
 public:
  PORTABLE_FORCEINLINE_FUNCTION
  Real *operator[](int i) const {
    return nullptr;
  }
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
  EOSChoice get() {
    return mpark::get<EOSChoice>(eos_);
  }

  Variant GetOnDevice() {
    return mpark::visit([](auto &eos) { return eos_variant<EOSs...>(eos.GetOnDevice()); },
                        eos_);
  }

  // Place member functions here
  PORTABLE_INLINE_FUNCTION
  Real TemperatureFromDensityInternalEnergy(const Real rho, const Real sie,
                                            Real *lambda = nullptr) const {
    return mpark::visit(
        [&rho, &sie, &lambda](const auto &eos) {
          return eos.TemperatureFromDensityInternalEnergy(rho, sie, lambda);
        },
        eos_);
  }

  PORTABLE_INLINE_FUNCTION
  Real InternalEnergyFromDensityTemperature(const Real rho, const Real temperature,
                                            Real *lambda = nullptr) const {
    return mpark::visit(
        [&rho, &temperature, &lambda](const auto &eos) {
          return eos.InternalEnergyFromDensityTemperature(rho, temperature, lambda);
        },
        eos_);
  }

  PORTABLE_INLINE_FUNCTION
  Real PressureFromDensityTemperature(const Real rho, const Real temperature,
                                      Real *lambda = nullptr) const {
    return mpark::visit(
        [&rho, &temperature, &lambda](const auto &eos) {
          return eos.PressureFromDensityTemperature(rho, temperature, lambda);
        },
        eos_);
  }

  PORTABLE_INLINE_FUNCTION
  Real PressureFromDensityInternalEnergy(const Real rho, const Real sie,
                                         Real *lambda = nullptr) const {
    return mpark::visit(
        [&rho, &sie, &lambda](const auto &eos) {
          return eos.PressureFromDensityInternalEnergy(rho, sie, lambda);
        },
        eos_);
  }

  PORTABLE_INLINE_FUNCTION
  Real SpecificHeatFromDensityTemperature(const Real rho, const Real temperature,
                                          Real *lambda = nullptr) const {
    return mpark::visit(
        [&rho, &temperature, &lambda](const auto &eos) {
          return eos.SpecificHeatFromDensityTemperature(rho, temperature, lambda);
        },
        eos_);
  }

  PORTABLE_INLINE_FUNCTION
  Real SpecificHeatFromDensityInternalEnergy(const Real rho, const Real sie,
                                             Real *lambda = nullptr) const {
    return mpark::visit(
        [&rho, &sie, &lambda](const auto &eos) {
          return eos.SpecificHeatFromDensityInternalEnergy(rho, sie, lambda);
        },
        eos_);
  }

  PORTABLE_INLINE_FUNCTION
  Real BulkModulusFromDensityTemperature(const Real rho, const Real temperature,
                                         Real *lambda = nullptr) const {
    return mpark::visit(
        [&rho, &temperature, &lambda](const auto &eos) {
          return eos.BulkModulusFromDensityTemperature(rho, temperature, lambda);
        },
        eos_);
  }

  PORTABLE_INLINE_FUNCTION
  Real BulkModulusFromDensityInternalEnergy(const Real rho, const Real sie,
                                            Real *lambda = nullptr) const {
    return mpark::visit(
        [&rho, &sie, &lambda](const auto &eos) {
          return eos.BulkModulusFromDensityInternalEnergy(rho, sie, lambda);
        },
        eos_);
  }

  PORTABLE_INLINE_FUNCTION
  Real GruneisenParamFromDensityTemperature(const Real rho, const Real temperature,
                                            Real *lambda = nullptr) const {
    return mpark::visit(
        [&rho, &temperature, &lambda](const auto &eos) {
          return eos.GruneisenParamFromDensityTemperature(rho, temperature, lambda);
        },
        eos_);
  }

  PORTABLE_INLINE_FUNCTION
  Real GruneisenParamFromDensityInternalEnergy(const Real rho, const Real sie,
                                               Real *lambda = nullptr) const {
    return mpark::visit(
        [&rho, &sie, &lambda](const auto &eos) {
          return eos.GruneisenParamFromDensityInternalEnergy(rho, sie, lambda);
        },
        eos_);
  }

  PORTABLE_INLINE_FUNCTION
  void FillEos(Real &rho, Real &temp, Real &energy, Real &press, Real &cv, Real &bmod,
               const unsigned long output, Real *lambda = nullptr) const {
    return mpark::visit(
        [&rho, &temp, &energy, &press, &cv, &bmod, &output, &lambda](const auto &eos) {
          return eos.FillEos(rho, temp, energy, press, cv, bmod, output, lambda);
        },
        eos_);
  }

  PORTABLE_INLINE_FUNCTION
  void ReferenceDensityTemperature(Real &rho, Real &T, Real *lambda = nullptr) const {
    return mpark::visit(
        [&rho, &T, &lambda](const auto &eos) {
          return eos.ReferenceDensityTemperature(rho, T, lambda);
        },
        eos_);
  }

  PORTABLE_INLINE_FUNCTION
  void ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                              Real &bmod, Real &dpde, Real &dvdt,
                              Real *lambda = nullptr) const {
    return mpark::visit(
        [&rho, &temp, &sie, &press, &cv, &bmod, &dpde, &dvdt, &lambda](const auto &eos) {
          return eos.ValuesAtReferenceState(rho, temp, sie, press, cv, bmod, dpde, dvdt,
                                            lambda);
        },
        eos_);
  }

  PORTABLE_INLINE_FUNCTION
  void PTofRE(Real &rho, Real &sie, Real *lambda, Real &press, Real &temp,
              Real &dpdr, Real &dpde, Real &dtdr, Real &dtde) const {
    return mpark::visit(
        [&rho, &sie, &lambda, &press, &temp, &dpdr, &dpde, &dtdr,
         &dtde](const auto &eos) {
          return eos.PTofRE(rho, sie, lambda, press, temp, dpdr, dpde, dtdr,
                            dtde);
        },
        eos_);
  }

  PORTABLE_INLINE_FUNCTION
  void DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                            Real *lambda, Real &rho, Real &sie) const {
    return mpark::visit(
        [&press, &temp, lambda, &rho, &sie](const auto &eos) {
          return eos.DensityEnergyFromPressureTemperature(press, temp, lambda, rho, sie);
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
  template<typename RealIndexer, typename ConstRealIndexer>
  inline
  void TemperatureFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&sies,
                                            RealIndexer &&temperatures,
                                            const int num) const {
    NullIndexer lambdas; // Returns null pointer for every index
    return TemperatureFromDensityInternalEnergy(
        &rhos, &sies, &temperatures, num, &lambdas);
  }

  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void TemperatureFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&sies,
                                            RealIndexer &&temperatures,
                                            const int num,
                                            LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &sies, &num, &temperatures, &num, &lambdas](const auto &eos) {
          return eos.TemperatureFromDensityInternalEnergy(
              rhos, sies, temperatures, num, lambdas);
        },
        eos_);
  }

  template<typename RealIndexer, typename ConstRealIndexer>
  inline
  void InternalEnergyFromDensityTemperature(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&temperatures,
                                            RealIndexer &&sies,
                                            const int num) const {
    NullIndexer lambdas; // Returns null pointer for every index
    return InternalEnergyFromDensityTemperature(
        &rhos, &temperatures, &sies, num, &lambdas);
  }

  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void InternalEnergyFromDensityTemperature(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&temperatures,
                                            RealIndexer &&sies,
                                            const int num,
                                            LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &temperatures, &sies, &num, &lambdas](const auto &eos) {
          return eos.InternalEnergyFromDensityTemperature(
              rhos, temperatures, sies, num, lambdas);
        },
        eos_);
  }

  template<typename RealIndexer, typename ConstRealIndexer>
  inline
  void PressureFromDensityTemperature(ConstRealIndexer &&rhos,
                                      ConstRealIndexer &&temperatures,
                                      RealIndexer &&pressures,
                                      const int num) const {
    NullIndexer lambdas; // Returns null pointer for every index
    return PressureFromDensityTemperature(
        &rhos, &temperatures, &pressures, num, &lambdas);
  }

  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void PressureFromDensityTemperature(ConstRealIndexer &&rhos,
                                      ConstRealIndexer &&temperatures,
                                      RealIndexer &&pressures,
                                      const int num,
                                      LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &temperatures, &pressures, &num, &lambdas](const auto &eos) {
          return eos.PressureFromDensityTemperature(
              rhos, temperatures, pressures, num, lambdas);
        },
        eos_);
  }

  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void PressureFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                         ConstRealIndexer &&sies,
                                         RealIndexer &&pressures,
                                         const int num,
                                         LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &sies, &pressures, &num, &lambdas](const auto &eos) {
          return eos.PressureFromDensityInternalEnergy(
              rhos, sies, pressures, num, lambdas);
        },
        eos_);
//    return mpark::visit(
//        [rhos, sies, pressures, num, lambdas](const auto &eos) {
//          return eos.PressureFromDensityInternalEnergy(
//              rhos, sies, pressures, num, lambdas);
//        },
//        eos_);
  }

  template<typename RealIndexer, typename ConstRealIndexer>
  inline
  void PressureFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                         ConstRealIndexer &&sies,
                                         RealIndexer &&pressures,
                                         const int num) const {
    NullIndexer lambdas{}; // Returns null pointer for every index
    return this->PressureFromDensityInternalEnergy(
        rhos, sies, pressures, num, lambdas);
  }

  template<typename RealIndexer, typename ConstRealIndexer>
  inline
  void SpecificHeatFromDensityTemperature(ConstRealIndexer &&rhos,
                                          ConstRealIndexer &&temperatures,
                                          RealIndexer &&cvs,
                                          const int num) const {
    NullIndexer lambdas; // Returns null pointer for every index
    return SpecificHeatFromDensityTemperature(
        &rhos, &temperatures, &cvs, num, &lambdas);
  }

  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void SpecificHeatFromDensityTemperature(ConstRealIndexer &&rhos,
                                          ConstRealIndexer &&temperatures,
                                          RealIndexer &&cvs,
                                          const int num,
                                          LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &temperatures, &cvs, &num, &lambdas](const auto &eos) {
          return eos.SpecificHeatFromDensityTemperature(
              rhos, temperatures, cvs, num, lambdas);
        },
        eos_);
  }

  template<typename RealIndexer, typename ConstRealIndexer>
  inline
  void SpecificHeatFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                             ConstRealIndexer &&sies,
                                             RealIndexer &&cvs,
                                             const int num) const {
    NullIndexer lambdas; // Returns null pointer for every index
    return SpecificHeatFromDensityInternalEnergy(
        &rhos, &sies, &cvs, num, &lambdas);
  }

  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void SpecificHeatFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                             ConstRealIndexer &&sies,
                                             RealIndexer &&cvs,
                                             const int num,
                                             LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &sies, &cvs, &num, &lambdas](const auto &eos) {
          return eos.SpecificHeatFromDensityInternalEnergy(
              rhos, sies, cvs, num, lambdas);
        },
        eos_);
  }

  template<typename RealIndexer, typename ConstRealIndexer>
  inline
  void BulkModulusFromDensityTemperature(ConstRealIndexer &&rhos,
                                         ConstRealIndexer &&temperatures,
                                         RealIndexer &&bmods,
                                         const int num) const {
    NullIndexer lambdas; // Returns null pointer for every index
    return BulkModulusFromDensityTemperature(
        &rhos, &temperatures, &bmods, num, &lambdas);
  }

  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void BulkModulusFromDensityTemperature(ConstRealIndexer &&rhos,
                                         ConstRealIndexer &&temperatures,
                                         RealIndexer &&bmods,
                                         const int num,
                                         LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &temperatures, &bmods, &num, &lambdas](const auto &eos) {
          return eos.BulkModulusFromDensityTemperature(
              rhos, temperatures, bmods, num, lambdas);
        },
        eos_);
  }

  template<typename RealIndexer, typename ConstRealIndexer>
  inline
  void BulkModulusFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&sies,
                                            RealIndexer &&bmods,
                                            const int num) const {
    NullIndexer lambdas; // Returns null pointer for every index
    return BulkModulusFromDensityInternalEnergy(
        &rhos, &sies, &bmods, num, &lambdas);
  }

  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void BulkModulusFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&sies,
                                            RealIndexer &&bmods,
                                            const int num,
                                            LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &sies, &bmods, &num, &lambdas](const auto &eos) {
          return eos.BulkModulusFromDensityInternalEnergy(
              rhos, sies, bmods, num, lambdas);
        },
        eos_);
  }

  template<typename RealIndexer, typename ConstRealIndexer>
  inline
  void GruneisenParamFromDensityTemperature(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&temperatures,
                                            RealIndexer &&gm1s,
                                            const int num) const {
    NullIndexer lambdas; // Returns null pointer for every index
    return GruneisenParamFromDensityTemperature(
        &rhos, &temperatures, &gm1s, num, &lambdas);
  }

  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void GruneisenParamFromDensityTemperature(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&temperatures,
                                            RealIndexer &&gm1s,
                                            const int num,
                                            LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &temperatures, &gm1s, &num, &lambdas](const auto &eos) {
          return eos.GruneisenParamFromDensityTemperature(
              rhos, temperatures, gm1s, num, lambdas);
        },
        eos_);
  }

  template<typename RealIndexer, typename ConstRealIndexer>
  inline
  void GruneisenParamFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                               ConstRealIndexer &&sies,
                                               RealIndexer &&gm1s,
                                               const int num) const {
    NullIndexer lambdas; // Returns null pointer for every index
    return GruneisenParamFromDensityInternalEnergy(
        &rhos, &sies, &gm1s, num, &lambdas);
  }

  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void GruneisenParamFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                               ConstRealIndexer &&sies,
                                               RealIndexer &&gm1s,
                                               const int num,
                                               LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &sies, &gm1s, &lambdas, &num](const auto &eos) {
          return eos.GruneisenParamFromDensityInternalEnergy(
              rhos, sies, gm1s, num, lambdas);
        },
        eos_);
  }

  template<typename RealIndexer>
  inline
  void FillEos(RealIndexer &&rhos, RealIndexer &&temps, RealIndexer &&energies,
               RealIndexer &&presses, RealIndexer &&cvs, RealIndexer &&bmods,
               const int num, const unsigned long output) const {
    NullIndexer lambdas; // Returns null pointer for every index
    return FillEos(&rhos, &temps, &energies, &presses, &cvs, &bmods,
                       num, output, &lambdas);
  }

  template<typename RealIndexer, typename LambdaIndexer>
  inline
  void FillEos(RealIndexer &&rhos, RealIndexer &&temps, RealIndexer &&energies,
               RealIndexer &&presses, RealIndexer &&cvs, RealIndexer &&bmods,
               const int num, const unsigned long output,
               LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &temps, &energies, &presses, &cvs, &bmods, &num, &output,
         &lambdas](const auto &eos) {
          return eos.FillEos(rhos, temps, energies, presses, cvs, bmods,
                                   num, output, lambdas);
        },
        eos_);
  }

  template<typename RealIndexer>
  inline
  void PTofRE(RealIndexer &&rhos, RealIndexer &&sies,
              RealIndexer &&presses, RealIndexer &&temps,
              RealIndexer &&dpdrs, RealIndexer &&dpdes,
              RealIndexer &&dtdrs, RealIndexer &&dtdes,
              const int num) const {
    NullIndexer lambdas; // Returns null pointer for every index
    return PTofRE(&rhos, &sies, &presses, &temps, &dpdrs, &dpdes, &dtdrs,
                  &dtdes, num, &lambdas);
  }

  template<typename RealIndexer, typename LambdaIndexer>
  inline
  void PTofRE(RealIndexer &&rhos, RealIndexer &&sies,
              RealIndexer &&presses, RealIndexer &&temps,
              RealIndexer &&dpdrs, RealIndexer &&dpdes,
              RealIndexer &&dtdrs, RealIndexer &&dtdes,
              const int num, LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &sies, &presses, &temps, &dpdrs, &dpdes, &dtdrs, &dtdes,
         &num, &lambdas](const auto &eos) {
          return eos.PTofRE(rhos, sies, presses, temps, dpdrs, dpdes, dtdrs,
                            dtdes, num, lambdas);
        },
        eos_);
  }

  PORTABLE_INLINE_FUNCTION
  unsigned long PreferredInput() const noexcept {
    return mpark::visit([](const auto &eos) { return eos.PreferredInput(); }, eos_);
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

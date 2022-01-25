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
#include <algorithm>

using Real = double;

namespace singularity {

template <typename... Ts>
using eos_variant = mpark::variant<Ts...>;

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
  void PTofRE(Real &rho, Real &sie, Real *lambda, Real &press, Real &temp, Real &dpdr,
              Real &dpde, Real &dtdr, Real &dtde) const {
    return mpark::visit(
        [&rho, &sie, &lambda, &press, &temp, &dpdr, &dpde, &dtdr,
         &dtde](const auto &eos) {
          press = eos.PressureFromDensityInternalEnergy(rho, sie, lambda);
          temp = eos.TemperatureFromDensityInternalEnergy(rho, sie, lambda);
          const Real drho = rho * 1.0e-6;
          const Real de = sie * 1.0e-6;
          const Real Pr = eos.PressureFromDensityInternalEnergy(rho + drho, sie, lambda);
          const Real Pe = eos.PressureFromDensityInternalEnergy(rho, sie + de, lambda);
          const Real Tr =
              eos.TemperatureFromDensityInternalEnergy(rho + drho, sie, lambda);
          const Real Te = eos.TemperatureFromDensityInternalEnergy(rho, sie + de, lambda);
          dpdr = (Pr - press) / drho;
          dpde = (Pe - press) / de;
          dtdr = (Tr - temp) / drho;
          dtde = (Te - temp) /
                 de; // Would it be better to skip the calculation of Te and return 1/cv?
          return;
        },
        eos_);
  }

  PORTABLE_INLINE_FUNCTION
  void DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                            Real *lambda, Real &rho, Real &sie) const {
    return mpark::visit(
        [&press, &temp, &lambda, &rho, &sie](const auto &eos) {
          return eos.DensityEnergyFromPressureTemperature(press, temp, lambda, rho, sie);
        },
        eos_);
  }

  /*
  Vector versions of the member functions are only available for CPU
  
  RealIndexer must have an operator[](int) that returns a Real. e.g., Real*
  ConstRealIndexer is as RealIndexer, but assumed const type.
  LambdaIndexer must have an operator[](int) that returns a Real*. e.g., Real**
*/
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void TemperatureFromDensityInternalEnergyVector(ConstRealIndexer &&rhos,
                                                  ConstRealIndexer &&sies,
                                                  const int num,
                                                  RealIndexer &&temperatures,
                                                  LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &sies, &num, &temperatures, &lambdas, &num](const auto &eos) {
          return eos.TemperatureFromDensityInternalEnergyVector(
              rhos, sies, num, temperatures, lambdas);
        },
        eos_);
  }

  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void InternalEnergyFromDensityTemperatureVector(ConstRealIndexer &&rhos,
                                                  ConstRealIndexer &&temperatures,
                                                  const int num,
                                                  RealIndexer &&sies,
                                                  LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &temperatures, &num, &sies, &lambdas](const auto &eos) {
          return eos.InternalEnergyFromDensityTemperatureVector(
              rhos, temperatures, sies, lambdas);
        },
        eos_);
  }

  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void PressureFromDensityTemperatureVector(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&temperatures,
                                            const int num,
                                            RealIndexer &&pressures,
                                            LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &temperatures, &num, &pressures, &lambdas](const auto &eos) {
          return eos.PressureFromDensityTemperatureVector(
              rhos, temperatures, pressures, lambdas);
        },
        eos_);
  }

  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void PressureFromDensityInternalEnergyVector(ConstRealIndexer &&rhos,
                                               ConstRealIndexer &&sies,
                                               const int num,
                                               RealIndexer &&pressures,
                                               LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &sies, &num, &pressures, &lambdas](const auto &eos) {
          return eos.PressureFromDensityInternalEnergyVector(
              rhos, sies, pressures, lambdas);
        },
        eos_);
  }

  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void SpecificHeatFromDensityTemperatureVector(ConstRealIndexer &&rhos,
                                                ConstRealIndexer &&temperatures,
                                                const int num,
                                                RealIndexer &&cvs,
                                                LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &temperatures, &num, &cvs, &lambdas](const auto &eos) {
          return eos.SpecificHeatFromDensityTemperatureVector(
              rhos, temperatures, cvs, lambdas);
        },
        eos_);
  }

  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void SpecificHeatFromDensityInternalEnergyVector(ConstRealIndexer &&rhos,
                                                   ConstRealIndexer &&sies,
                                                   const int num,
                                                   RealIndexer &&cvs,
                                                   LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &sies, &num, &cvs, &lambdas](const auto &eos) {
          return eos.SpecificHeatFromDensityInternalEnergyVector(
              rhos, sies, cvs, lambdas);
        },
        eos_);
  }

  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void BulkModulusFromDensityTemperatureVector(ConstRealIndexer &&rhos,
                                               ConstRealIndexer &&temperatures,
                                               const int num,
                                               RealIndexer &&bmods,
                                               LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &temperatures, &num, &bmods, &lambdas](const auto &eos) {
          return eos.BulkModulusFromDensityTemperatureVector(
              rhos, temperatures, bmods, lambdas);
        },
        eos_);
  }

  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void BulkModulusFromDensityInternalEnergyVector(ConstRealIndexer &&rhos,
                                                  ConstRealIndexer &&sies,
                                                  const int num,
                                                  RealIndexer &&bmods,
                                                  LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &sies, &num, &bmods, &lambdas](const auto &eos) {
          return eos.BulkModulusFromDensityInternalEnergyVector(
              rhos, sies, bmods, lambdas);
        },
        eos_);
  }

  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void GruneisenParamFromDensityTemperatureVector(ConstRealIndexer &&rhos,
                                                  ConstRealIndexer &&temperatures,
                                                  const int num,
                                                  RealIndexer &&gm1s,
                                                  LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &temperatures, &num, &gm1s, &lambdas](const auto &eos) {
          return eos.GruneisenParamFromDensityTemperatureVector(
              rhos, temperatures, gm1s, lambdas);
        },
        eos_);
  }

  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void GruneisenParamFromDensityInternalEnergyVector(ConstRealIndexer &&rhos,
                                                     ConstRealIndexer &&sies,
                                                     const int num,
                                                     RealIndexer &&gm1s,
                                                     LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &sies, &num, &gm1s, &lambdas](const auto &eos) {
          return eos.GruneisenParamFromDensityInternalEnergyVector(
              rhos, sies, gm1s, lambdas);
        },
        eos_);
  }

  template<typename RealIndexer, typename LambdaIndexer>
  inline
  void FillEosVector(RealIndexer &&rhos, RealIndexer &&temps, RealIndexer &&energies,
                     RealIndexer &&presses, RealIndexer &&cvs, RealIndexer &&bmods,
                     const int num, const unsigned long output,
                     LambdaIndexer &&lambdas) const {
    return mpark::visit(
        [&rhos, &temps, &energies, &presses, &cvs, &bmods, &num, &output,
         &lambdas](const auto &eos) {
          return eos.FillEosVector(rhos, temps, energies, presses, cvs, bmods,
                                   num, output, lambdas);
        },
        eos_);
  }

  template<typename RealIndexer, typename LambdaIndexer>
  inline
  void PTofREVector(RealIndexer &&rhos, RealIndexer &&sies, const int num,
                    LambdaIndexer &&lambdas, RealIndexer &&presses,
                    RealIndexer &&temps, RealIndexer &&dpdrs,
                    RealIndexer &&dpdes, RealIndexer &&dtdrs,
                    RealIndexer &&dtdes) const {
    return mpark::visit(
        [&rhos, &sies, &num, &lambdas, &presses, &temps, &dpdrs, &dpdes, &dtdrs,
         &dtdes](const auto &eos) {

          // Lookup at density and energy
          eos.PressureFromDensityInternalEnergyVector(rhos, sies, num, presses,
                                                      lambdas);
          eos.TemperatureFromDensityInternalEnergyVector(rhos, sies, num, temps,
                                                         lambdas);
          // Peturbation factor
          Real factor = 1. + 1.0e-06;

          // Perturb densities and do lookups
          for (int i = 0; i < num ; i ++) {
            rhos[i] *= factor;
          }
          eos.PressureFromDensityInternalEnergyVector(rhos, sies, dpdrs, lambdas);
          eos.TemperatureFromDensityInternalEnergyVector(rhos, sies, dtdrs, lambdas);
          
          // Reset densities, perturb energies, and do lookups
          for (int i = 0; i < num ; i ++) {
            sies[i] *= factor;
            rhos[i] /= factor;
          }
          eos.PressureFromDensityInternalEnergyVector(rhos, sies, dpdrs, lambdas);
          eos.TemperatureFromDensityInternalEnergyVector(rhos, sies, dtdrs, lambdas);

          // Reset the energies to their original values
          for (int i = 0; i < num ; i ++) {
            sies[i] /= factor;
          }

          // Calculate the derivatives
          for(int i = 0; i < num ; i ++){
            dpdrs[i] = (dpdrs[i] - presses[i]) / (rhos[i] * (1. - factor));
            dpdes[i] = (dpdes[i] - presses[i]) / (sies[i] * (1. - factor));
            dtdrs[i] = (dtdrs[i] - temps[i]) / (rhos[i] * (1. - factor));
            dtdes[i] = (dtdes[i] - temps[i]) / (sies[i] * (1. - factor));
          }

          return;
        },
        eos_);
  }

  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void DensityEnergyFromPressureTemperatureVector(ConstRealIndexer &&presses,
                                                  ConstRealIndexer &&temps,
                                                  LambdaIndexer &&lambdas,
                                                  RealIndexer &rhos,
                                                  RealIndexer &sies) const {
    return mpark::visit(
        [&presses, &temps, &lambdas, &rhos, &sies](const auto &eos) {
          return eos.DensityEnergyFromPressureTemperatureVector(
            presses, temps, lambdas, rhos, sies);
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

//------------------------------------------------------------------------------
// © 2021-2026. Triad National Security, LLC. All rights reserved. This
// program was produced under U.S. Government contract 89233218CNA000001
// for Los Alamos National Laboratory (LANL), which is operated by Triad
// National Security, LLC for the U.S. Department of Energy/National
// Nuclear Security Administration. All rights in the program are
// reserved by Triad National Security, LLC, and the U.S. Department of
// Energy/National Nuclear Security Administration. The Government is
// granted for itself and others acting on its behalf a nonexclusive,
// paid-up, irrevocable worldwide license in this material to reproduce,
// prepare derivative works, distribute copies to the public, perform
// publicly and display publicly, and to permit others to do so.
//------------------------------------------------------------------------------

#ifndef _SINGULARITY_EOS_EOS_FREE_ENERGY_SCALE_
#define _SINGULARITY_EOS_EOS_FREE_ENERGY_SCALE_

#include <cmath>
#include <cstdio>
#include <string>
#include <utility>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_errors.hpp>
#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/eos_error.hpp>
#include <singularity-eos/eos/eos_base.hpp>
#include <singularity-eos/eos/modifiers/modifier_vector_macros.hpp>

namespace singularity {

using namespace eos_base;

// A thermodynamically consistent, multiplicative scaling of the specific
// Helmholtz free energy:
//
//   f_new(rho, T) = amplitude * f_base(rho, T).
//
// Consequently, pressure, internal energy, entropy, specific heat, and bulk
// modulus all scale by amplitude. The Gruneisen parameter is unchanged. The
// transformation is valid only while amplitude is positive and fixed with
// respect to thermodynamic derivatives.
template <typename T>
class FreeEnergyAmplitudeScaling : public EosBase<FreeEnergyAmplitudeScaling<T>> {
 public:
  SG_ADD_BASE_CLASS_USINGS(FreeEnergyAmplitudeScaling<T>);
  using BaseType = T;

  static std::string EosType() {
    return std::string("FreeEnergyAmplitudeScaling<") + T::EosType() + std::string(">");
  }

  static std::string EosPyType() {
    return std::string("FreeEnergyAmplitudeScaling") + T::EosPyType();
  }

  // Move semantics ensures dynamic memory comes along for the ride.
  FreeEnergyAmplitudeScaling(T &&t, const Real amplitude)
      : t_(std::forward<T>(t)), amplitude_(amplitude), inv_amplitude_(1. / amplitude) {
    CheckParams();
  }
  FreeEnergyAmplitudeScaling() = default;

  PORTABLE_INLINE_FUNCTION void CheckParams() const {
    PORTABLE_ALWAYS_REQUIRE(!std::isnan(amplitude_),
                            "Free-energy amplitude must be well defined.");
    PORTABLE_ALWAYS_REQUIRE(!std::isnan(inv_amplitude_),
                            "Inverse free-energy amplitude must be well defined.");
    PORTABLE_ALWAYS_REQUIRE(amplitude_ > 0, "Free-energy amplitude must be positive.");
    PORTABLE_ALWAYS_REQUIRE(inv_amplitude_ > 0,
                            "Inverse free-energy amplitude must be positive.");
    t_.CheckParams();
  }

  auto GetOnDevice() {
    return FreeEnergyAmplitudeScaling<T>(t_.GetOnDevice(), amplitude_);
  }
  inline void Finalize() { t_.Finalize(); }

  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    return t_.TemperatureFromDensityInternalEnergy(rho, inv_amplitude_ * sie, lambda);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temperature, Indexer_t &&lambda = nullptr) const {
    return amplitude_ * t_.InternalEnergyFromDensityTemperature(rho, temperature, lambda);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real PressureFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    return amplitude_ *
           t_.PressureFromDensityInternalEnergy(rho, inv_amplitude_ * sie, lambda);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real
  MinInternalEnergyFromDensity(const Real rho, Indexer_t &&lambda = nullptr) const {
    return amplitude_ * t_.MinInternalEnergyFromDensity(rho, lambda);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real EntropyFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    return amplitude_ *
           t_.EntropyFromDensityInternalEnergy(rho, inv_amplitude_ * sie, lambda);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    return amplitude_ *
           t_.SpecificHeatFromDensityInternalEnergy(rho, inv_amplitude_ * sie, lambda);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    return amplitude_ *
           t_.BulkModulusFromDensityInternalEnergy(rho, inv_amplitude_ * sie, lambda);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    return t_.GruneisenParamFromDensityInternalEnergy(rho, inv_amplitude_ * sie, lambda);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real PressureFromDensityTemperature(
      const Real rho, const Real temperature, Indexer_t &&lambda = nullptr) const {
    return amplitude_ * t_.PressureFromDensityTemperature(rho, temperature, lambda);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real EntropyFromDensityTemperature(
      const Real rho, const Real temperature, Indexer_t &&lambda = nullptr) const {
    return amplitude_ * t_.EntropyFromDensityTemperature(rho, temperature, lambda);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real SpecificHeatFromDensityTemperature(
      const Real rho, const Real temperature, Indexer_t &&lambda = nullptr) const {
    return amplitude_ * t_.SpecificHeatFromDensityTemperature(rho, temperature, lambda);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real BulkModulusFromDensityTemperature(
      const Real rho, const Real temperature, Indexer_t &&lambda = nullptr) const {
    return amplitude_ * t_.BulkModulusFromDensityTemperature(rho, temperature, lambda);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temperature, Indexer_t &&lambda = nullptr) const {
    return t_.GruneisenParamFromDensityTemperature(rho, temperature, lambda);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION void
  InternalEnergyFromDensityPressure(const Real rho, const Real press, Real &sie,
                                    Indexer_t &&lambda = nullptr) const {
    sie *= inv_amplitude_;
    t_.InternalEnergyFromDensityPressure(rho, inv_amplitude_ * press, sie, lambda);
    sie *= amplitude_;
  }

  template <typename Lambda_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  PTDerivativesFromPreferred(const Real rho, const Real sie, const Real press,
                             const Real temp, Lambda_t &&lambda, Real &dedP_T,
                             Real &drdP_T, Real &dedT_P, Real &drdT_P) const {
    t_.PTDerivativesFromPreferred(rho, inv_amplitude_ * sie, inv_amplitude_ * press, temp,
                                  lambda, dedP_T, drdP_T, dedT_P, drdT_P);
    // de_new/dP_new at fixed T is unchanged. Density derivatives with
    // respect to pressure scale inversely, while energy derivatives with
    // respect to temperature scale with the amplitude.
    drdP_T *= inv_amplitude_;
    dedT_P *= amplitude_;
  }

  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION void FillEos(Real &rho, Real &temp, Real &energy, Real &press,
                                 Real &cv, Real &bmod, const unsigned long output,
                                 Indexer_t &&lambda = nullptr) const {
    Real senergy =
        (output & thermalqs::specific_internal_energy) ? 0 : inv_amplitude_ * energy;
    Real spress = (output & thermalqs::pressure) ? 0 : inv_amplitude_ * press;
    Real scv = 0;
    Real sbmod = 0;
    t_.FillEos(rho, temp, senergy, spress, scv, sbmod, output, lambda);
    if (output & thermalqs::specific_internal_energy) energy = amplitude_ * senergy;
    if (output & thermalqs::pressure) press = amplitude_ * spress;
    if (output & thermalqs::specific_heat) cv = amplitude_ * scv;
    if (output & thermalqs::bulk_modulus) bmod = amplitude_ * sbmod;
  }

  template <typename Indexer_t = Real *>
  PORTABLE_FUNCTION void ValuesAtReferenceState(Real &rho, Real &temp, Real &sie,
                                                Real &press, Real &cv, Real &bmod,
                                                Real &dpde, Real &dvdt,
                                                Indexer_t &&lambda = nullptr) const {
    t_.ValuesAtReferenceState(rho, temp, sie, press, cv, bmod, dpde, dvdt, lambda);
    sie *= amplitude_;
    press *= amplitude_;
    cv *= amplitude_;
    bmod *= amplitude_;
    // dpde and dvdt are invariant under an amplitude scaling.
  }

  // Vector implementations
  SG_MODIFIER_FORWARD_2IN_1OUT(TemperatureFromDensityInternalEnergy,
                               transform.y.apply(inv_amplitude_);)
  SG_MODIFIER_FORWARD_2IN_1OUT(PressureFromDensityTemperature,
                               transform.f.apply(amplitude_);)
  SG_MODIFIER_FORWARD_2IN_1OUT(PressureFromDensityInternalEnergy,
                               transform.y.apply(inv_amplitude_);
                               transform.f.apply(amplitude_);)
  SG_MODIFIER_FORWARD_1IN_1OUT(MinInternalEnergyFromDensity,
                               transform.f.apply(amplitude_);)
  SG_MODIFIER_FORWARD_2IN_1OUT(SpecificHeatFromDensityTemperature,
                               transform.f.apply(amplitude_);)
  SG_MODIFIER_FORWARD_2IN_1OUT(SpecificHeatFromDensityInternalEnergy,
                               transform.y.apply(inv_amplitude_);
                               transform.f.apply(amplitude_);)
  SG_MODIFIER_FORWARD_2IN_1OUT(BulkModulusFromDensityTemperature,
                               transform.f.apply(amplitude_);)
  SG_MODIFIER_FORWARD_2IN_1OUT(BulkModulusFromDensityInternalEnergy,
                               transform.y.apply(inv_amplitude_);
                               transform.f.apply(amplitude_);)
  SG_MODIFIER_FORWARD_2IN_1OUT(GruneisenParamFromDensityTemperature,
                               transform.x.apply(1.);)
  SG_MODIFIER_FORWARD_2IN_1OUT(GruneisenParamFromDensityInternalEnergy,
                               transform.y.apply(inv_amplitude_);)
  SG_MODIFIER_FORWARD_2IN_1OUT(InternalEnergyFromDensityTemperature,
                               transform.f.apply(amplitude_);)
  SG_MODIFIER_FORWARD_2IN_1OUT(EntropyFromDensityTemperature,
                               transform.f.apply(amplitude_);)
  SG_MODIFIER_FORWARD_2IN_1OUT(EntropyFromDensityInternalEnergy,
                               transform.y.apply(inv_amplitude_);
                               transform.f.apply(amplitude_);)

  constexpr static inline int nlambda() noexcept { return T::nlambda(); }

  template <typename Indexable>
  static inline constexpr bool NeedsLambda() {
    return T::template NeedsLambda<Indexable>();
  }

  static constexpr unsigned long PreferredInput() { return T::PreferredInput(); }

  static inline unsigned long scratch_size(std::string method, unsigned int nelements) {
    return T::scratch_size(method, nelements);
  }

  static inline unsigned long max_scratch_size(unsigned int nelements) {
    return T::max_scratch_size(nelements);
  }

  PORTABLE_FUNCTION void PrintParams() const {
    t_.PrintParams();
    printf("free_energy_amplitude = %f\n", amplitude_);
  }

  template <typename Indexer_t>
  PORTABLE_FUNCTION void
  DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                       Indexer_t &&lambda, Real &rho, Real &sie) const {
    t_.DensityEnergyFromPressureTemperature(inv_amplitude_ * press, temp, lambda, rho,
                                            sie);
    sie *= amplitude_;
  }

  PORTABLE_FORCEINLINE_FUNCTION Real MinimumDensity() const {
    return t_.MinimumDensity();
  }

  PORTABLE_FORCEINLINE_FUNCTION Real MinimumTemperature() const {
    return t_.MinimumTemperature();
  }

  PORTABLE_FORCEINLINE_FUNCTION Real MaximumDensity() const {
    return t_.MaximumDensity();
  }

  PORTABLE_FORCEINLINE_FUNCTION Real MinimumPressure() const {
    return amplitude_ * t_.MinimumPressure();
  }

  PORTABLE_FORCEINLINE_FUNCTION Real MaximumPressureAtTemperature(const Real temp) const {
    return amplitude_ * t_.MaximumPressureAtTemperature(temp);
  }

  PORTABLE_FORCEINLINE_FUNCTION Real RhoPmin(const Real temp) const {
    return t_.RhoPmin(temp);
  }

  PORTABLE_FORCEINLINE_FUNCTION Real MinimumInternalEnergy() const {
    return amplitude_ * t_.MinimumInternalEnergy();
  }

  PORTABLE_FORCEINLINE_FUNCTION Real MaximumInternalEnergy() const {
    return amplitude_ * t_.MaximumInternalEnergy();
  }

  PORTABLE_INLINE_FUNCTION Real MeanAtomicMass() const { return t_.MeanAtomicMass(); }

  PORTABLE_INLINE_FUNCTION Real MeanAtomicNumber() const { return t_.MeanAtomicNumber(); }

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

 private:
  T t_;
  double amplitude_;
  double inv_amplitude_;
};

} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_FREE_ENERGY_SCALE_

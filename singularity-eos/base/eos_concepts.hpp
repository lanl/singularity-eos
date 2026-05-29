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

#ifndef _SINGULARITY_EOS_BASE_EOS_CONCEPTS_HPP_
#define _SINGULARITY_EOS_BASE_EOS_CONCEPTS_HPP_

// This file was generated in part with the assistance of generative AI

namespace eos_concepts {

template <typename EOS>
inline constexpr bool dependent_false_v = false;

// C++17 implementation using detection idiom
#if __cplusplus < 202002L

// Detection helper
template <typename...>
using void_t = void;

// Detect MinimumDensity
template <typename EOS, typename = void>
struct has_MinimumDensity : std::false_type {};

template <typename EOS>
struct has_MinimumDensity<EOS, 
                          void_t<decltype(std::declval<EOS>().MinimumDensity()))>> : std::true_type {};

// Detect MaximumDensity
template <typename EOS, typename = void>
struct has_MaximumDensity : std::false_type {};

template <typename EOS>
struct has_MaximumDensity<EOS, 
                          void_t<decltype(std::declval<EOS>().MaximumDensity()))>> : std::true_type {};             
// Detect MinimumTemperature
template <typename EOS, typename = void>
struct has_MinimumTemperature : std::false_type {};

template <typename EOS>
struct has_MinimumTemperature<EOS, 
                          void_t<decltype(std::declval<EOS>().MinimumTemperature()))>> : std::true_type {};   

// Detect PressureFromDensityTemperature
template <typename EOS, typename = void>
struct has_P_rho_T : std::false_type {};

template <typename EOS>
struct has_P_rho_T<EOS,
                   void_t<decltype(std::declval<EOS>().PressureFromDensityTemperature(
                       std::declval<Real>(), std::declval<Real>()))>> : std::true_type {};

// Detect PressureFromDensityInternalEnergy
template <typename EOS, typename = void>
struct has_P_rho_sie : std::false_type {};

template <typename EOS>
struct has_P_rho_sie<
    EOS, void_t<decltype(std::declval<EOS>().PressureFromDensityInternalEnergy(
             std::declval<Real>(), std::declval<Real>()))>> : std::true_type {};

// Detect TemperatureFromDensityInternalEnergy
template <typename EOS, typename = void>
struct has_T_rho_sie : std::false_type {};

template <typename EOS>
struct has_T_rho_sie<
    EOS, void_t<decltype(std::declval<EOS>().TemperatureFromDensityInternalEnergy(
             std::declval<Real>(), std::declval<Real>()))>> : std::true_type {};

// Detect InternalEnergyFromDensityTemperature
template <typename EOS, typename = void>
struct has_sie_rho_T : std::false_type {};

template <typename EOS>
struct has_sie_rho_T<
    EOS, void_t<decltype(std::declval<EOS>().InternalEnergyFromDensityTemperature(
             std::declval<Real>(), std::declval<Real>()))>> : std::true_type {};

// Detect BulkModulusFromDensityTemperature
template <typename EOS, typename = void>
struct has_bmod_rho_T : std::false_type {};

template <typename EOS>
struct has_bmod_rho_T<
    EOS, void_t<decltype(std::declval<EOS>().BulkModulusFromDensityTemperature(
             std::declval<Real>(), std::declval<Real>()))>> : std::true_type {};

// Detect BulkModulusFromDensityInternalEnergy
template <typename EOS, typename = void>
struct has_bmod_rho_sie : std::false_type {};

template <typename EOS>
struct has_bmod_rho_sie<
    EOS, void_t<decltype(std::declval<EOS>().BulkModulusFromDensityInternalEnergy(
             std::declval<Real>(), std::declval<Real>()))>> : std::true_type {};

// Detect GruneisenParamFromDensityTemperature
template <typename EOS, typename = void>
struct has_gamma_rho_T : std::false_type {};

template <typename EOS>
struct has_gamma_rho_T<
    EOS, void_t<decltype(std::declval<EOS>().GruneisenParamFromDensityTemperature(
             std::declval<Real>(), std::declval<Real>()))>> : std::true_type {};

// Detect GruneisenParamFromDensityInternalEnergy
template <typename EOS, typename = void>
struct has_gamma_rho_sie : std::false_type {};

template <typename EOS>
struct has_gamma_rho_sie<
    EOS, void_t<decltype(std::declval<EOS>().GruneisenParamFromDensityInternalEnergy(
             std::declval<Real>(), std::declval<Real>()))>> : std::true_type {};

// Detect SpecificHeatFromDensityTemperature
template <typename EOS, typename = void>
struct has_cv_rho_T : std::false_type {};

template <typename EOS>
struct has_cv_rho_T<
    EOS, void_t<decltype(std::declval<EOS>().SpecificHeatFromDensityTemperature(
             std::declval<Real>(), std::declval<Real>()))>> : std::true_type {};

// Detect SpecificHeatFromDensityInternalEnergy
template <typename EOS, typename = void>
struct has_cv_rho_sie : std::false_type {};

template <typename EOS>
struct has_cv_rho_sie<
    EOS, void_t<decltype(std::declval<EOS>().SpecificHeatFromDensityInternalEnergy(
             std::declval<Real>(), std::declval<Real>()))>> : std::true_type {};

// Detect MeanAtomicMass() method
template <typename EOS, typename = void>
struct has_abar : std::false_type {};

template <typename EOS>
struct has_abar<EOS, void_t<decltype(std::declval<EOS>().MeanAtomicMass())>>
    : std::true_type {};

// Detect MeanAtomicNumber() method
template <typename EOS, typename = void>
struct has_zbar : std::false_type {};

template <typename EOS>
struct has_zbar<EOS, void_t<decltype(std::declval<EOS>().MeanAtomicNumber())>>
    : std::true_type {};

// Convenience constexpr bools

template <typename EOS>
inline constexpr bool has_MinimumDensity_v = has_MinimumDensity<EOS>.value;

template <typename EOS>
inline constexpr bool has_MaximumDensity_v = has_MaximumDensity<EOS>.value;

template <typename EOS>
inline constexpr bool has_MinimumTemperature_v = has_MinimumTemperature<EOS>.value;

template <typename EOS>
inline constexpr bool has_P_rho_T_v = has_P_rho_T<EOS>::value;

template <typename EOS>
inline constexpr bool has_P_rho_sie_v = has_P_rho_sie<EOS>::value;

template <typename EOS>
inline constexpr bool has_T_rho_sie_v = has_T_rho_sie<EOS>::value;

template <typename EOS>
inline constexpr bool has_sie_rho_T_v = has_sie_rho_T<EOS>::value;

template <typename EOS>
inline constexpr bool has_bmod_rho_T_v = has_bmod_rho_T<EOS>::value;

template <typename EOS>
inline constexpr bool has_bmod_rho_sie_v = has_bmod_rho_sie<EOS>::value;

template <typename EOS>
inline constexpr bool has_gamma_rho_T_v = has_gamma_rho_T<EOS>::value;

template <typename EOS>
inline constexpr bool has_gamma_rho_sie_v = has_gamma_rho_sie<EOS>::value;

template <typename EOS>
inline constexpr bool has_cv_rho_T_v = has_cv_rho_T<EOS>::value;

template <typename EOS>
inline constexpr bool has_cv_rho_sie_v = has_cv_rho_sie<EOS>::value;

template <typename EOS>
inline constexpr bool has_abar_v = has_abar<EOS>::value;

template <typename EOS>
inline constexpr bool has_zbar_v = has_zbar<EOS>::value;

#else
// C++20 implementation using concepts (for future migration)

template <typename EOS>
concept has_MinimumDensity = requires(EOS eos) {
  { eos.MinimumDensity } -> std::same_as<Real>;
};

template <typename EOS>
concept has_MaximumDensity = requires(EOS eos) {
  { eos.MaximumDensity } -> std::same_as<Real>;
};

template <typename EOS>
concept has_MinimumTemperature = requires(EOS eos) {
  { eos.MinimumTemperature } -> std::same_as<Real>;
};

template <typename EOS>
concept has_P_rho_T = requires(EOS eos, Real rho, Real T) {
  { eos.PressureFromDensityTemperature(rho, T) } -> std::same_as<Real>;
};

template <typename EOS>
concept has_P_rho_sie = requires(EOS eos, Real rho, Real sie) {
  { eos.PressureFromDensityInternalEnergy(rho, sie) } -> std::same_as<Real>;
};

template <typename EOS>
concept has_T_rho_sie = requires(EOS eos, Real rho, Real sie) {
  { eos.TemperatureFromDensityInternalEnergy(rho, sie) } -> std::same_as<Real>;
};

template <typename EOS>
concept has_sie_rho_T = requires(EOS eos, Real rho, Real T) {
  { eos.InternalEnergyFromDensityTemperature(rho, T) } -> std::same_as<Real>;
};

template <typename EOS>
concept has_bmod_rho_T = requires(EOS eos, Real rho, Real T) {
  { eos.BulkModulusFromDensityTemperature(rho, T) } -> std::same_as<Real>;
};

template <typename EOS>
concept has_bmod_rho_sie = requires(EOS eos, Real rho, Real sie) {
  { eos.BulkModulusFromDensityInternalEnergy(rho, sie) } -> std::same_as<Real>;
};

template <typename EOS>
concept has_gamma_rho_T = requires(EOS eos, Real rho, Real T) {
  { eos.GruneisenParamFromDensityTemperature(rho, T) } -> std::same_as<Real>;
};

template <typename EOS>
concept has_gamma_rho_sie = requires(EOS eos, Real rho, Real sie) {
  { eos.GruneisenParamFromDensityInternalEnergy(rho, sie) } -> std::same_as<Real>;
};

template <typename EOS>
concept has_cv_rho_T = requires(EOS eos, Real rho, Real T) {
  { eos.SpecificHeatFromDensityTemperature(rho, T) } -> std::same_as<Real>;
};

template <typename EOS>
concept has_cv_rho_sie = requires(EOS eos, Real rho, Real sie) {
  { eos.SpecificHeatFromDensityInternalEnergy(rho, sie) } -> std::same_as<Real>;
};

template <typename EOS>
concept has_abar = requires(EOS eos) {
  { eos.MeanAtomicMass() } -> std::same_as<Real>;
};

template <typename EOS>
concept has_zbar = requires(EOS eos) {
  { eos.MeanAtomicNumber() } -> std::same_as<Real>;
};

// For compatibility with C++17 code, provide _v helpers
template <typename EOS>
inline constexpr bool has_MinimumDensity_v = has_MinimumDensity<EOS>;

template <typename EOS>
inline constexpr bool has_MinimumTemperature_v = has_MinimumTemperature<EOS>;

template <typename EOS>
inline constexpr bool has_MaximumDensity_v = has_MaximumDensity<EOS>;

template <typename EOS>
inline constexpr bool has_P_rho_T_v = has_P_rho_T<EOS>;

template <typename EOS>
inline constexpr bool has_P_rho_sie_v = has_P_rho_sie<EOS>;

template <typename EOS>
inline constexpr bool has_T_rho_sie_v = has_T_rho_sie<EOS>;

template <typename EOS>
inline constexpr bool has_sie_rho_T_v = has_sie_rho_T<EOS>;

template <typename EOS>
inline constexpr bool has_bmod_rho_T_v = has_bmod_rho_T<EOS>;

template <typename EOS>
inline constexpr bool has_bmod_rho_sie_v = has_bmod_rho_sie<EOS>;

template <typename EOS>
inline constexpr bool has_cv_rho_T_v = has_cv_rho_T<EOS>;

template <typename EOS>
inline constexpr bool has_cv_rho_sie_v = has_cv_rho_sie<EOS>;

template <typename EOS>
inline constexpr bool has_gamma_rho_T_v = has_gamma_rho_T<EOS>;

template <typename EOS>
inline constexpr bool has_gamma_rho_sie_v = has_gamma_rho_sie<EOS>;

template <typename EOS>
inline constexpr bool has_abar_v = has_abar<EOS>;

template <typename EOS>
inline constexpr bool has_zbar_v = has_zbar<EOS>;

#endif

} // namespace eos_concepts

#endif // SINGULARITY_EOS_BASE_EOS_CONCEPTS_HPP_
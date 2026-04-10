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

// This file was created in part with generative AI

#ifndef _SINGULARITY_EOS_EOS_MODIFIERS_MODIFIER_VECTOR_MACROS_
#define _SINGULARITY_EOS_EOS_MODIFIERS_MODIFIER_VECTOR_MACROS_

#include <utility>

#include <singularity-eos/eos/eos_base.hpp>

/* TODO(JMM):
   These macros shrink the amount of boiler plate overloads necessary
   for modifiers.
 */
#define SG_MODIFIER_DEVICE_WRAP_2IN_1OUT(NAME)                                           \
  template <typename LambdaIndexer>                                                      \
  inline void NAME(const Real *in1, const Real *in2, Real *out, Real *scratch,           \
                   const int num, LambdaIndexer &&lambdas,                               \
                   Transform &&transform = Transform()) const {                          \
    NAME(PortsOfCall::Exec::Device(), in1, in2, out, scratch, num,                       \
         std::forward<LambdaIndexer>(lambdas), std::forward<Transform>(transform));      \
  }

#define SG_MODIFIER_DEVICE_WRAP_1IN_1OUT(NAME)                                           \
  template <typename LambdaIndexer>                                                      \
  inline void NAME(const Real *in1, Real *out, Real *scratch, const int num,             \
                   LambdaIndexer &&lambdas, Transform &&transform = Transform()) const { \
    NAME(PortsOfCall::Exec::Device(), in1, out, scratch, num,                            \
         std::forward<LambdaIndexer>(lambdas), std::forward<Transform>(transform));      \
  }

#define SG_MODIFIER_DEVICE_WRAP_ALL()                                                    \
  SG_MODIFIER_DEVICE_WRAP_2IN_1OUT(TemperatureFromDensityInternalEnergy)                 \
  SG_MODIFIER_DEVICE_WRAP_2IN_1OUT(PressureFromDensityTemperature)                       \
  SG_MODIFIER_DEVICE_WRAP_2IN_1OUT(PressureFromDensityInternalEnergy)                    \
  SG_MODIFIER_DEVICE_WRAP_1IN_1OUT(MinInternalEnergyFromDensity)                         \
  SG_MODIFIER_DEVICE_WRAP_2IN_1OUT(SpecificHeatFromDensityTemperature)                   \
  SG_MODIFIER_DEVICE_WRAP_2IN_1OUT(SpecificHeatFromDensityInternalEnergy)                \
  SG_MODIFIER_DEVICE_WRAP_2IN_1OUT(BulkModulusFromDensityTemperature)                    \
  SG_MODIFIER_DEVICE_WRAP_2IN_1OUT(BulkModulusFromDensityInternalEnergy)                 \
  SG_MODIFIER_DEVICE_WRAP_2IN_1OUT(GruneisenParamFromDensityTemperature)                 \
  SG_MODIFIER_DEVICE_WRAP_2IN_1OUT(GruneisenParamFromDensityInternalEnergy)              \
  SG_MODIFIER_DEVICE_WRAP_2IN_1OUT(InternalEnergyFromDensityTemperature)                 \
  SG_MODIFIER_DEVICE_WRAP_2IN_1OUT(EntropyFromDensityTemperature)                        \
  SG_MODIFIER_DEVICE_WRAP_2IN_1OUT(EntropyFromDensityInternalEnergy)

#define SG_MODIFIER_FORWARD_2IN_1OUT(NAME, PREPARE)                                      \
  template <typename Space, typename LambdaIndexer,                                      \
            typename EnableIfSpace =                                                     \
                std::enable_if_t<!variadic_utils::has_int_index_v<Space>>>               \
  inline void NAME(const Space &s, const Real *in1, const Real *in2, Real *out,          \
                   Real *scratch, const int num, LambdaIndexer &&lambdas,                \
                   Transform &&transform = Transform()) const {                          \
    PREPARE                                                                              \
    t_.NAME(s, in1, in2, out, scratch, num, std::forward<LambdaIndexer>(lambdas),        \
            std::forward<Transform>(transform));                                         \
  }                                                                                      \
  SG_MODIFIER_DEVICE_WRAP_2IN_1OUT(NAME)

#define SG_MODIFIER_FORWARD_1IN_1OUT(NAME, PREPARE)                                      \
  template <typename Space, typename LambdaIndexer,                                      \
            typename EnableIfSpace =                                                     \
                std::enable_if_t<!variadic_utils::has_int_index_v<Space>>>               \
  inline void NAME(const Space &s, const Real *in1, Real *out, Real *scratch,            \
                   const int num, LambdaIndexer &&lambdas,                               \
                   Transform &&transform = Transform()) const {                          \
    PREPARE                                                                              \
    t_.NAME(s, in1, out, scratch, num, std::forward<LambdaIndexer>(lambdas),             \
            std::forward<Transform>(transform));                                         \
  }                                                                                      \
  SG_MODIFIER_DEVICE_WRAP_1IN_1OUT(NAME)

#endif

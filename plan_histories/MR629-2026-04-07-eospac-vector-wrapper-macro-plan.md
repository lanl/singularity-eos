# MR629 EOSPAC Vector Wrapper Macro Plan

## Goal

Reduce the boilerplate in `singularity-eos/eos/eos_eospac.hpp` for the repeated EOSPAC vector wrapper declarations that forward `2-in/1-out` calls through `EosBase<EOSPAC>`.

## Scope

Only target the repeated wrapper pattern with:

- two indexed real inputs
- one indexed real output
- the no-scratch overload
- the scratch overload used for type-mismatch fallback

Do not introduce a macro for the warning strings themselves. Keep both warnings written verbatim inside the wrapper macro body.

Do not add a `1-in/1-out` macro. `MinInternalEnergyFromDensity` stays as-is.

Do not touch the specialized pointer-based EOSPAC vector implementations, `FillEos`, `DensityEnergyFromPressureTemperature`, or `ValuesAtReferenceState`.

## Proposed Macro

Add one EOSPAC-local macro near the vector wrapper declarations:

```cpp
#define SG_EOSPAC_VEC_2IN_1OUT(NAME, IN1, IN2, OUT)                                      \
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>     \
  inline void NAME(ConstRealIndexer &&IN1, ConstRealIndexer &&IN2, RealIndexer &&OUT,    \
                   const int num, LambdaIndexer &&lambdas) const {                       \
    PORTABLE_WARN("Not providing scratch memory will trigger scalar EOSPAC lookups");    \
    EosBase<EOSPAC>::NAME(std::forward<ConstRealIndexer>(IN1),                           \
                          std::forward<ConstRealIndexer>(IN2),                           \
                          std::forward<RealIndexer>(OUT), num,                           \
                          std::forward<LambdaIndexer>(lambdas));                         \
  }                                                                                      \
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer,     \
            typename = std::enable_if_t<!is_raw_pointer<RealIndexer, Real>::value>>      \
  inline void NAME(ConstRealIndexer &&IN1, ConstRealIndexer &&IN2, RealIndexer &&OUT,    \
                   Real * /*scratch*/, const int num, LambdaIndexer &&lambdas) const {   \
    PORTABLE_WARN("EOSPAC type mismatch will cause significant performance degradation"); \
    EosBase<EOSPAC>::NAME(std::forward<ConstRealIndexer>(IN1),                           \
                          std::forward<ConstRealIndexer>(IN2),                           \
                          std::forward<RealIndexer>(OUT), num,                           \
                          std::forward<LambdaIndexer>(lambdas));                         \
  }
```

## Planned Replacements

Use the macro for these wrapper declarations:

```cpp
SG_EOSPAC_VEC_2IN_1OUT(TemperatureFromDensityInternalEnergy, rhos, sies, temperatures)
SG_EOSPAC_VEC_2IN_1OUT(InternalEnergyFromDensityTemperature, rhos, temperatures, sies)
SG_EOSPAC_VEC_2IN_1OUT(PressureFromDensityTemperature, rhos, temperatures, pressures)
SG_EOSPAC_VEC_2IN_1OUT(PressureFromDensityInternalEnergy, rhos, sies, pressures)
SG_EOSPAC_VEC_2IN_1OUT(EntropyFromDensityTemperature, rhos, temperatures, entropies)
SG_EOSPAC_VEC_2IN_1OUT(EntropyFromDensityInternalEnergy, rhos, sies, entropies)
SG_EOSPAC_VEC_2IN_1OUT(SpecificHeatFromDensityTemperature, rhos, temperatures, cvs)
SG_EOSPAC_VEC_2IN_1OUT(SpecificHeatFromDensityInternalEnergy, rhos, sies, cvs)
SG_EOSPAC_VEC_2IN_1OUT(BulkModulusFromDensityTemperature, rhos, temperatures, bmods)
SG_EOSPAC_VEC_2IN_1OUT(BulkModulusFromDensityInternalEnergy, rhos, sies, bmods)
SG_EOSPAC_VEC_2IN_1OUT(GruneisenParamFromDensityTemperature, rhos, temperatures, gm1s)
SG_EOSPAC_VEC_2IN_1OUT(GruneisenParamFromDensityInternalEnergy, rhos, sies, gm1s)
```

## Rationale

This follows the same general style as the vector API boilerplate macros in `singularity-eos/eos/eos_base.hpp`, but keeps the implementation local to EOSPAC because the warning behavior is EOSPAC-specific.

The result should remove most of the duplicated declarations without obscuring the actual fallback behavior.

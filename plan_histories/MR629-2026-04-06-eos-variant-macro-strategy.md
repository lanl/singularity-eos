# Plan: Reduce `Variant` Boilerplate in `eos_variant.hpp`

## Goal

Reduce duplication in `singularity-eos/eos/eos_variant.hpp` by introducing a small set of helper macros for the repeated `PortsOfCall::visit(...)` wrappers, without changing the public API of `Variant`.

## Summary

A macro-based cleanup should work well in `eos_variant.hpp`.

The vector section is especially repetitive and appears to be a better fit for macros than the already-refactored `EosBase` vector wrappers. Most methods repeat the same four-overload bundle:

1. a no-lambda overload that constructs `NullIndexer`
2. a lambda-taking overload that dispatches through `PortsOfCall::visit`
3. a `scratch` + no-lambda overload that constructs `NullIndexer`
4. a `scratch` + lambda-taking overload that dispatches through `PortsOfCall::visit`

## Recommended Macro Families

### 1. Vector wrappers

Use a small family of macros for the repeated vector forwarding patterns:

- `SG_VARIANT_VEC_2IN_1OUT(NAME, IN1, IN2, OUT)`
- `SG_VARIANT_VEC_1IN_1OUT(NAME, IN1, OUT)`
- `SG_VARIANT_VEC_2IN_REFOUT(NAME, IN1, IN2, OUT)`

These would cover the methods that forward to the same method on the active EOS in the variant.

### 2. Optional scalar wrappers

If further cleanup is desired, a second macro family could reduce the scalar `PortsOfCall::visit(...)` wrappers:

- `SG_VARIANT_SCALAR_2IN_1OUT(NAME, A1, A2)`
- `SG_VARIANT_SCALAR_1IN_1OUT(NAME, A1)`
- `SG_VARIANT_SCALAR_2IN_REFOUT(NAME, A1, A2, OUT)`

This scalar cleanup is optional. The highest-value, lowest-risk target is the vector section.

## Best Targets For Vector Macros

These methods appear to fit the macro pattern well:

- `TemperatureFromDensityInternalEnergy`
- `InternalEnergyFromDensityTemperature`
- `PressureFromDensityTemperature`
- `PressureFromDensityInternalEnergy`
- `MinInternalEnergyFromDensity`
- `EntropyFromDensityTemperature`
- `EntropyFromDensityInternalEnergy`
- `GibbsFreeEnergyFromDensityTemperature`
- `GibbsFreeEnergyFromDensityInternalEnergy`
- `SpecificHeatFromDensityTemperature`
- `SpecificHeatFromDensityInternalEnergy`
- `BulkModulusFromDensityTemperature`
- `BulkModulusFromDensityInternalEnergy`
- `GruneisenParamFromDensityTemperature`
- `GruneisenParamFromDensityInternalEnergy`
- `InternalEnergyFromDensityPressure`

## Methods To Leave Handwritten

These are less regular and should likely remain explicit:

- `FillEos`
- `ReferenceDensityTemperature`
- `ValuesAtReferenceState`
- `DensityEnergyFromPressureTemperature`
- various introspection helpers
- serialization helpers

## Suggested Shape

The vector section could collapse to a short list of macro invocations such as:

```cpp
SG_VARIANT_VEC_2IN_1OUT(TemperatureFromDensityInternalEnergy, rhos, sies, temperatures)
SG_VARIANT_VEC_2IN_1OUT(InternalEnergyFromDensityTemperature, rhos, temperatures, sies)
SG_VARIANT_VEC_2IN_1OUT(PressureFromDensityTemperature, rhos, temperatures, pressures)
SG_VARIANT_VEC_2IN_1OUT(PressureFromDensityInternalEnergy, rhos, sies, pressures)
SG_VARIANT_VEC_1IN_1OUT(MinInternalEnergyFromDensity, rhos, sies)
SG_VARIANT_VEC_2IN_1OUT(EntropyFromDensityTemperature, rhos, temperatures, entropies)
SG_VARIANT_VEC_2IN_1OUT(EntropyFromDensityInternalEnergy, rhos, sies, entropies)
SG_VARIANT_VEC_2IN_1OUT(GibbsFreeEnergyFromDensityTemperature, rhos, temperatures, Gs)
SG_VARIANT_VEC_2IN_1OUT(GibbsFreeEnergyFromDensityInternalEnergy, rhos, sies, Gs)
SG_VARIANT_VEC_2IN_1OUT(SpecificHeatFromDensityTemperature, rhos, temperatures, cvs)
SG_VARIANT_VEC_2IN_1OUT(SpecificHeatFromDensityInternalEnergy, rhos, sies, cvs)
SG_VARIANT_VEC_2IN_1OUT(BulkModulusFromDensityTemperature, rhos, temperatures, bmods)
SG_VARIANT_VEC_2IN_1OUT(BulkModulusFromDensityInternalEnergy, rhos, sies, bmods)
SG_VARIANT_VEC_2IN_1OUT(GruneisenParamFromDensityTemperature, rhos, temperatures, gm1s)
SG_VARIANT_VEC_2IN_1OUT(GruneisenParamFromDensityInternalEnergy, rhos, sies, gm1s)
SG_VARIANT_VEC_2IN_REFOUT(InternalEnergyFromDensityPressure, rhos, Ps, sies)
```

## Important Constraint

The macros in `eos_variant.hpp` will be slightly more delicate than those in `eos_base.hpp`, because they use perfect forwarding through `PortsOfCall::visit` lambdas. That is still manageable, but it argues for a small, focused macro family rather than an overly generic one.

## Recommendation

Proceed with the vector macro refactor first.

It offers the best reduction in boilerplate with limited risk and preserves the existing overload set. If the result reads well, consider a second pass for the scalar `visit(...)` wrappers.

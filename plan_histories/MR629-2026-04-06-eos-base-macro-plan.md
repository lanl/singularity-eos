# Plan: Reduce `EosBase` Vector Wrapper Boilerplate with Macros

## Goal

Reduce duplication in `singularity-eos/eos/eos_base.hpp` for the CRTP-based vector API wrappers that call `portableFor`, without fundamentally changing the public API or overload set.

## Scope

Target the repeated vector wrapper patterns in `EosBase<CRTP>`:

- one-input / one-output wrappers
- two-input / one-output wrappers
- wrappers with a `scratch` forwarding overload
- wrappers with a raw-pointer forwarding overload
- the special vector wrapper that forwards to a scalar method with an output reference

Leave non-matching methods handwritten if the macro would become harder to read than the original code.

## Proposed Macro Shapes

### 1. `SG_EOS_VEC_2IN_1OUT`

Use for methods shaped like:

```cpp
out[i] = copy.Method(in1[i], in2[i], lambdas[i]);
```

Examples:

- `TemperatureFromDensityInternalEnergy`
- `InternalEnergyFromDensityTemperature`
- `PressureFromDensityTemperature`
- `PressureFromDensityInternalEnergy`
- `EntropyFromDensityTemperature`
- `EntropyFromDensityInternalEnergy`
- `SpecificHeatFromDensityTemperature`
- `SpecificHeatFromDensityInternalEnergy`
- `BulkModulusFromDensityTemperature`
- `BulkModulusFromDensityInternalEnergy`
- `GruneisenParamFromDensityTemperature`
- `GruneisenParamFromDensityInternalEnergy`
- `GibbsFreeEnergyFromDensityTemperature`
- `GibbsFreeEnergyFromDensityInternalEnergy`

Suggested parameters:

- method name
- first indexed input variable name
- second indexed input variable name
- output variable name
- copy declaration kind: `const CRTP &` or `CRTP`

### 2. `SG_EOS_VEC_1IN_1OUT`

Use for methods shaped like:

```cpp
out[i] = copy.Method(in1[i], lambdas[i]);
```

Primary example:

- `MinInternalEnergyFromDensity`

Suggested parameters:

- method name
- indexed input variable name
- output variable name
- copy declaration kind

### 3. `SG_EOS_VEC_2IN_REFOUT`

Use for methods shaped like:

```cpp
copy.Method(in1[i], in2[i], out[i], lambdas[i]);
```

Primary example:

- `InternalEnergyFromDensityPressure`

Suggested parameters:

- method name
- first indexed input variable name
- second indexed input variable name
- output variable name
- copy declaration kind

## Important Constraints

- Preserve the current overload names and signatures.
- Preserve the `scratch` overload behavior.
- Preserve the raw-pointer overload behavior with `Transform && = Transform()`.
- Preserve current `copy` semantics:
  some methods use `const CRTP &copy`, while others use `CRTP copy`.
- Keep `FillEos` handwritten unless a separate dedicated macro remains clearly readable.

## Suggested Implementation Order

1. Add the helper macros near the other EOS helper macros at the top of `eos_base.hpp`.
2. Replace the repetitive vector wrapper bodies with macro invocations.
3. Keep `FillEos` and any non-conforming methods handwritten.
4. Build and run the existing test suite that covers EOS vector calls.
5. Check compiler diagnostics carefully, since macro errors can be noisy in templated code.

## Expected Benefits

- substantially less boilerplate in `EosBase`
- easier addition of new vector wrapper methods
- lower risk of copy-paste inconsistencies between overloads
- no fundamental change to the API seen by derived EOS classes

## Main Risks

- macro-generated template errors may be harder to read
- overly generic macros can become less maintainable than the duplicated code
- accidental mismatch between current `const CRTP &copy` and `CRTP copy` behavior

## Recommendation

Proceed with a small macro family rather than one fully generic macro. Three focused macros should remove most duplication while keeping the generated code readable and close to the current implementation.

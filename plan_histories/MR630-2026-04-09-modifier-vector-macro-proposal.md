# MR630 Modifier Vector Macro Proposal

## Goal

Reduce the line count in `singularity-eos/eos/modifiers/*.hpp` for the raw-pointer
vector API, especially the overloads that only forward to
`PortsOfCall::Exec::Device()`.

## Current Situation

The modified vector API in the modifier headers now has two overload families for each
raw-pointer vector method:

1. A `Space` overload:

```cpp
template <typename Space, typename LambdaIndexer,
          typename EnableIfSpace =
              std::enable_if_t<!variadic_utils::has_int_index_v<Space>>>
inline void NAME(const Space &s, ...);
```

2. A default-device wrapper:

```cpp
template <typename LambdaIndexer>
inline void NAME(...) const {
  NAME(PortsOfCall::Exec::Device(), ...);
}
```

This repeats across:

- `scaled_eos.hpp`
- `shifted_eos.hpp`
- `eos_unitsystem.hpp`
- `ramps_eos.hpp`
- `floored_energy.hpp`

## Important Constraint

The default-device wrappers cannot simply be deleted while preserving the current API.
The inherited wrappers in `eos/eos_base.hpp` do not automatically route into the
modifier-specific raw-pointer `Space` overloads. If we want to keep the existing
call surface, the realistic near-term win is to generate these wrappers, not remove the
behavior.

## Recommended Macro Scheme

Introduce a small modifier-only macro header, for example:

- `singularity-eos/eos/modifiers/modifier_vector_macros.hpp`

That header would define the default-device wrappers once.

### 2-in / 1-out wrapper

```cpp
#define SG_MODIFIER_DEVICE_WRAP_2IN_1OUT(NAME)                                           \
  template <typename LambdaIndexer>                                                      \
  inline void NAME(const Real *in1, const Real *in2, Real *out, Real *scratch,           \
                   const int num, LambdaIndexer &&lambdas,                               \
                   Transform &&transform = Transform()) const {                          \
    NAME(PortsOfCall::Exec::Device(), in1, in2, out, scratch, num,                       \
         std::forward<LambdaIndexer>(lambdas), std::forward<Transform>(transform));      \
  }
```

### 1-in / 1-out wrapper

```cpp
#define SG_MODIFIER_DEVICE_WRAP_1IN_1OUT(NAME)                                           \
  template <typename LambdaIndexer>                                                      \
  inline void NAME(const Real *in1, Real *out, Real *scratch, const int num,             \
                   LambdaIndexer &&lambdas,                                              \
                   Transform &&transform = Transform()) const {                          \
    NAME(PortsOfCall::Exec::Device(), in1, out, scratch, num,                            \
         std::forward<LambdaIndexer>(lambdas), std::forward<Transform>(transform));      \
  }
```

### Aggregate wrapper list

```cpp
#define SG_MODIFIER_DEVICE_WRAP_ALL()                                                    \
  SG_MODIFIER_DEVICE_WRAP_2IN_1OUT(TemperatureFromDensityInternalEnergy)                  \
  SG_MODIFIER_DEVICE_WRAP_2IN_1OUT(PressureFromDensityTemperature)                        \
  SG_MODIFIER_DEVICE_WRAP_2IN_1OUT(PressureFromDensityInternalEnergy)                     \
  SG_MODIFIER_DEVICE_WRAP_1IN_1OUT(MinInternalEnergyFromDensity)                          \
  SG_MODIFIER_DEVICE_WRAP_2IN_1OUT(SpecificHeatFromDensityTemperature)                    \
  SG_MODIFIER_DEVICE_WRAP_2IN_1OUT(SpecificHeatFromDensityInternalEnergy)                 \
  SG_MODIFIER_DEVICE_WRAP_2IN_1OUT(BulkModulusFromDensityTemperature)                     \
  SG_MODIFIER_DEVICE_WRAP_2IN_1OUT(BulkModulusFromDensityInternalEnergy)                  \
  SG_MODIFIER_DEVICE_WRAP_2IN_1OUT(GruneisenParamFromDensityTemperature)                  \
  SG_MODIFIER_DEVICE_WRAP_2IN_1OUT(GruneisenParamFromDensityInternalEnergy)               \
  SG_MODIFIER_DEVICE_WRAP_2IN_1OUT(InternalEnergyFromDensityTemperature)                  \
  SG_MODIFIER_DEVICE_WRAP_2IN_1OUT(EntropyFromDensityTemperature)                         \
  SG_MODIFIER_DEVICE_WRAP_2IN_1OUT(EntropyFromDensityInternalEnergy)
```

## How This Would Be Used

Each modifier would keep only its `Space` overload implementations. At the end of the
vector API block, it would add:

```cpp
SG_MODIFIER_DEVICE_WRAP_ALL()
```

That removes most of the duplicate default-device wrapper code while preserving the
current behavior.

## Second Stage for Simple Forwarders

For `scaled_eos.hpp` and `eos_unitsystem.hpp`, the `Space` overload bodies are also
highly repetitive:

- update `transform`
- call `t_.NAME(s, ...)`

Those two files could use a second macro family such as:

```cpp
#define SG_MODIFIER_FORWARD_2IN_1OUT(NAME, PREPARE)                                      \
  template <typename Space, typename LambdaIndexer,                                      \
            typename EnableIfSpace =                                                     \
                std::enable_if_t<!variadic_utils::has_int_index_v<Space>>>               \
  inline void NAME(const Space &s, const Real *in1, const Real *in2, Real *out,          \
                   Real *scratch, const int num, LambdaIndexer &&lambdas,                \
                   Transform &&transform = Transform()) const {                          \
    PREPARE                                                                              \
    t_.NAME(s, in1, in2, out, scratch, num,                                              \
            std::forward<LambdaIndexer>(lambdas), std::forward<Transform>(transform));   \
  }                                                                                      \
  SG_MODIFIER_DEVICE_WRAP_2IN_1OUT(NAME)
```

Example usage:

```cpp
SG_MODIFIER_FORWARD_2IN_1OUT(
    TemperatureFromDensityInternalEnergy,
    transform.x.apply(scale_);
    transform.y.apply(inv_scale_);)
```

## Where Not to Overuse Macros

I would not force the body macros into:

- `shifted_eos.hpp`
- `floored_energy.hpp`
- `ramps_eos.hpp`

Those files contain real per-method logic:

- scratch-buffer preprocessing and postprocessing
- local `portableFor` loops
- ramp-pressure overrides

For those files, the better tradeoff is:

- keep explicit `Space` overload bodies
- macro-generate only the default-device wrappers

## Alternative If True Elimination Is Desired

If the actual goal is to eliminate per-file default-device overloads entirely, use a CRTP
mixin instead of macros. Each modifier would define only the `Space` overloads, and the
mixin would provide the no-space wrappers once.

That is conceptually cleaner than making the macro scheme more aggressive, but it is a
larger refactor than the macro-only approach above.

## Recommendation

Preferred implementation order:

1. Add a shared modifier macro header that generates the default-device wrappers.
2. Convert all five modified modifier headers to use that wrapper macro set.
3. Optionally add the forwarding-body macros for `ScaledEOS` and `EOSUnitSystem`.
4. Leave `ShiftedEOS`, `FlooredEnergy`, and `BilinearRampEOS` with explicit `Space`
   overload bodies.

This gives most of the line-count reduction without hiding the custom vector logic that
still needs to be readable.

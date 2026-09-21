# Plan: Expose Internal-Energy Bounds Through Modifiers (starting with `UnitSystem`)

**MR**: #XXX (rename this file once the MR number is assigned)
**Date**: 2026-09-19
**Status**: ✅ Phase 1 implemented — 📋 Phase 2 proposed, not implemented

## Motivation

A host code uses

```cpp
using baseEOS = singularity::SpinerEOSDependsRhoSie;
using EOS     = singularity::UnitSystem<baseEOS>;
```

and needs the table's energy bounds to seed a root-finder / bounds array. Today the only
way to get them is to strip the modifier:

```cpp
auto unmodified   = eosh.GetUnmodifiedObject();
sie_bounds(1)             = unmodified.sieMin();
sie_bounds(nbounds+1)     = unmodified.sieMax();
```

which returns values in the *base* (cgs) unit system, not in the unit system the host
code is actually working in. Every use site therefore has to remember to multiply by the
energy unit by hand, which is exactly the kind of conversion `UnitSystem` exists to
absorb.

`UnitSystem` already converts the other introspection quantities
(`eos_unitsystem.hpp:275-293`: `MinimumDensity`, `MinimumTemperature`, `MaximumDensity`,
`MinimumPressure`, `MaximumPressureAtTemperature`, `RhoPmin`). Energy bounds are simply
missing from that set — for `UnitSystem` and for the introspection API in general.

## Current Situation

### Where bounds introspection lives

- `EosBase` supplies permissive defaults for every model:
  `eos_base.hpp:535` (`MinimumDensity() -> 0`), `:537` (`MinimumTemperature() -> 0`),
  `:546` (`MaximumDensity() -> 1e100`), `:551` (`MinimumPressure() -> 0`),
  `:554` (`MaximumPressureAtTemperature() -> 1e100`), `:558` (`RhoPmin() -> 0`).
  There is **no** energy-bound member in this set.
- Modifiers either forward verbatim via `SG_ADD_MODIFIER_INTROSPECTION_METHODS`
  (`eos_base.hpp:133`; used by `shifted_eos.hpp:429`, `floored_energy.hpp:413`,
  `relativistic_eos.hpp:210`, `zsplit_eos.hpp:283`) or hand-write transformed versions
  (`scaled_eos.hpp:253-269`, `eos_unitsystem.hpp:275-293`, `ramps_eos.hpp:253-269`).
- `Variant` re-exposes each one through a `PortsOfCall::visit` (`eos_variant.hpp:475-503`).

### Where energy bounds exist today

- `SpinerEOSDependsRhoSie::sieMin()/sieMax()` — `eos_spiner_rho_sie.hpp:265-270`
- `StellarCollapse::sieMin()/sieMax()` — `eos_stellar_collapse.hpp:225-226`
- EOSPAC knows them (`eospac_wrapper.cpp:70-71` fills `metadata.sieMin/sieMax`) but
  `EOSPAC` itself only stores `rho_min_`/`temp_min_` (`eos_eospac.hpp:1212-1213`).
- `SpinerEOSDependsRhoT` and `Helmholtz` expose density/temperature bounds only
  (`eos_spiner_rho_temp.hpp:237-255`, `eos_helmholtz.hpp:269-273`).

These are ad-hoc, model-specific accessors — not part of the `EosBase` contract — which
is why no modifier and no variant can see them.

## Two Scopes

### Phase 1 — `UnitSystem` pass-through (small, unblocks the host code)

Add to `eos_unitsystem.hpp`, adjacent to the existing introspection block (~line 293):

```cpp
  PORTABLE_FORCEINLINE_FUNCTION Real sieMin() const { return inv_sie_unit_ * t_.sieMin(); }
  PORTABLE_FORCEINLINE_FUNCTION Real sieMax() const { return inv_sie_unit_ * t_.sieMax(); }
```

`inv_sie_unit_` is the correct factor: the file's convention is "multiply by a unit to
convert to cgs", which is why `MinInternalEnergyFromDensity` (`eos_unitsystem.hpp:132`)
returns `inv_sie_unit_ * S`.

**Why this does not require touching any other EOS.** `UnitSystem<T>` is a class
template, so member function *bodies* are only instantiated when called. Implicit
instantiation of `UnitSystem<IdealGas>` — including as a variant alternative — instantiates
member *declarations* only. The code above therefore compiles for every `T` in the
variant, and only fails if someone actually calls `sieMin()` on a `UnitSystem<T>` whose
`T` lacks it. That failure is a clear compile error at the call site, not a silent wrong
answer.

Optional, same cost, for symmetry with the base tables: `rhoMin()`, `rhoMax()`,
`TMin()`, `TMax()` pass-throughs (scaled by `inv_rho_unit_` / `inv_temp_unit_`). These
partly duplicate `MinimumDensity`/`MaximumDensity`/`MinimumTemperature`, so include them
only if the host code wants the base-table naming to survive the modifier.

**Explicit limitation of Phase 1**: this works on the concrete
`UnitSystem<SpinerEOSDependsRhoSie>` type only. It does *not* work through
`singularity::EOS`:

- `Variant::GetUnmodifiedObject()` returns a `Variant`, not a concrete EOS
  (`eos_variant.hpp:797`), so `.sieMin()` on the result does not compile.
- Inside a `visit` / `EvaluateHost` / `EvaluateDevice` lambda the body must compile for
  *every* alternative, so a bare `eos.sieMin()` fails on `IdealGas`. A host-side
  `if constexpr` + detection trait works but pushes a sentinel-value convention into
  downstream code — i.e. a worse version of Phase 2.

#### Phase 1 Implementation Record (2026-09-19)

Implemented as described above, energy bounds only; the optional
`rhoMin`/`rhoMax`/`TMin`/`TMax` pass-throughs were **not** added, since
`MinimumDensity`/`MaximumDensity`/`MinimumTemperature` already cover them through the
modifier.

Files changed:

- `singularity-eos/eos/modifiers/eos_unitsystem.hpp`: `sieMin()`/`sieMax()` after the
  `RhoPmin` introspection block, with a comment recording the lazy-instantiation
  requirement.
- `test/test_eos_modifiers.cpp`: `sieMin()/sieMax()` added to the `BoundedGas` helper;
  new `THEN` block in the "Modifiers propagate introspection bounds correctly" scenario
  asserting the `sie_unit` conversion and the round trip back to base units; comment on
  the `UnitSystem<IdealGas>` scenario noting that it guards usability for a `T` without
  the accessors.
- `doc/sphinx/src/modifiers.rst`: documented both methods in the unit system section,
  including the compile-time-error caveat and the replaced `GetUnmodifiedObject` idiom.
- `CHANGELOG.md`: entry under `## Current develop` → `### Added` (PR number is a
  placeholder).

Verification:

- `test/test_eos_modifiers.cpp` built and run standalone against the system Catch2:
  all tests pass, 26 assertions in 3 test cases.
- Separate syntax-only compile of `UnitSystem<SpinerEOSDependsRhoSie>::sieMin()/sieMax()`
  (the host-code pattern that motivated this) and of `UnitSystem<IdealGas>` with only
  `MinimumDensity()` called, confirming the lazy-instantiation property in both
  directions.
- `clang-format` (v21.1.4) applied; no reflow of surrounding code.

Not verified: the full CMake test suite. The pre-existing `build/` directory cannot
reconfigure because its spack-installed `ports-of-call`/`spiner` prefixes under `/tmp`
have been removed, so the standalone compiles above were used instead.

### Phase 2 — first-class energy bounds in the introspection API (optional, larger)

Only needed if the bounds must be reachable from a runtime-polymorphic
`singularity::EOS`.

Proposed names: `MinimumInternalEnergy()` / `MaximumInternalEnergy()`. Deliberately not
`MinimumEnergy`, to avoid confusion with the existing density-dependent cold curve
`MinInternalEnergyFromDensity`, which is a different quantity (a curve, not a table
bound). Naming is a review decision — flag it early.

1. **`EosBase` defaults** (`eos_base.hpp`, next to `:535-551`):
   `MinimumInternalEnergy() -> -1e100`, `MaximumInternalEnergy() -> 1e100`. Note the
   default minimum must be *negative*: energies are legitimately negative for cold curves
   and shifted EOS, so `0` is not a safe floor. Follow the existing `MaximumDensity`
   comment (`eos_base.hpp:539-545`) on big-finite-number vs. infinity.
2. **`SG_ADD_MODIFIER_INTROSPECTION_METHODS`** (`eos_base.hpp:133`): add verbatim
   forwarding of both methods. This covers `RelativisticEOS` (energy passes through
   untouched, `relativistic_eos.hpp:72-73`) and `RampEOS` (modifies pressure only)
   correctly for free.
3. **Modifiers needing a real transform** — these must *not* use the plain macro:
   - `UnitSystem`: `inv_sie_unit_ * t_.M*InternalEnergy()`.
   - `ScaledEOS`: modified energy is `scale_ * base` (`scaled_eos.hpp:73,80`), so
     `scale_ * t_.M*InternalEnergy()`.
   - `ShiftedEOS`: modified energy is `base + shift_` (`shifted_eos.hpp:75,85`). It
     currently uses the plain macro (`shifted_eos.hpp:429`), so it needs an explicit
     override; forwarding unchanged would be wrong by exactly `shift_`.
4. **Modifiers with a documented caveat rather than a transform**:
   - `ZSplit` scales energy by a lambda-dependent factor (`zsplit_eos.hpp:65,88`), and the
     introspection API takes no lambda. Forward unchanged and document that the bound is
     the un-split bound.
   - `FlooredEnergy` clamps energy to the cold curve per-density; the global bound is
     unchanged, so forwarding is correct.
5. **Concrete models**: override in `SpinerEOSDependsRhoSie` and `StellarCollapse`
   (trivial — delegate to existing `sieMin()/sieMax()`). `EOSPAC` can store the values
   already computed in `eospac_wrapper.cpp:70-71`. `SpinerEOSDependsRhoT` and `Helmholtz`
   are open questions (see below) and can keep the base defaults initially.
6. **`Variant`**: add two `visit`-based accessors alongside `eos_variant.hpp:475-503`.
7. **Python bindings**: add to the variant bindings; `SpinerEOSDependsRhoSie` and
   `StellarCollapse` already expose `sieMin`/`sieMax` (`python/module.cpp:143-144,157-158`).

## Sign Caveat (applies to both phases, worth noting in review)

`ScaledEOS::CheckParams` only requires `|scale| > 0` (`scaled_eos.hpp:59-65`), so a
negative scale is legal. Multiplying a min bound by a negative number turns it into a max
bound. `ScaledEOS::MinimumDensity` (`scaled_eos.hpp:253`) already has this latent issue,
so Phase 2 should either `min`/`max`-swap when `scale_ < 0` or explicitly document that
negative scales are unsupported for bounds introspection. `UnitSystem` is safe — its
`CheckParams` requires strictly positive units (`eos_unitsystem.hpp:100-105`).

## Testing

- `test/test_eos_modifiers.cpp` already has a `BoundedGas` helper with hand-written bounds
  (lines 64-82) and unit-system bound assertions (lines 328-375). Extend `BoundedGas` with
  energy bounds and add:
  - `UnitSystem`: energy bounds divided by `sie_unit`.
  - Phase 2 only: `ShiftedEOS` shifts by `shift_`, `ScaledEOS` scales by `scale_`,
    unmodified analytic EOS returns the permissive defaults, and the same values come back
    through `singularity::EOS`.
- Phase 1 needs at least one test that instantiates `UnitSystem<T>` for a `T` *without*
  `sieMin` (e.g. `IdealGas`) and exercises the rest of its API, to lock in the
  lazy-instantiation guarantee the design relies on.
- Tabulated coverage (`test/test_eos_tabulated.cpp`) for the Spiner round trip:
  `UnitSystem<SpinerEOSDependsRhoSie>(...).sieMin() * sie_unit == base.sieMin()`.

## Documentation

- `doc/sphinx/src/using-eos.rst`, "Methods Used for Mixed Cell Closures" (~lines
  1490-1535): document the new methods next to `MinimumDensity`/`MaximumDensity`, and
  extend the existing `warning` block about unbounded EOS to cover the energy defaults.
- `doc/sphinx/src/modifiers.rst`: note how each modifier transforms energy bounds.

## Bookkeeping (repo checklist)

- `CHANGELOG.md` under `## Current develop` → `### Added`, with the `[[PRxxx]]` link form.
- Update copyright year on each modified file; `eos_unitsystem.hpp` already carries the
  generative-AI notice (line 15), other touched files need one added if AI-assisted.
- `make format` after configuring.
- Rename this file to match the MR number.

## Open Questions

1. Ship Phase 1 alone, or Phase 1 + Phase 2 together? Phase 1 is ~6 lines and unblocks the
   host code immediately; Phase 2 is the API-consistent answer but touches ~10 files and
   needs review on naming and defaults.
2. `MinimumInternalEnergy` vs. `MinimumEnergy` vs. `SieMin` for the Phase 2 names.
3. `SpinerEOSDependsRhoT` energy bounds: derive from the tabulated `sie` field's min/max,
   or from `InternalEnergyFromDensityTemperature` at the table corners? The corner
   evaluation is only valid if `sie` is monotone in both arguments over the table.
4. Should `EOSPAC` store the `sieMin`/`sieMax` it already computes at load time?
5. Aside: `SG_ADD_MODIFIER_INTROSPECTION_METHODS(t)` (`eos_base.hpp:133`) names its
   parameter `t` but its body uses `t_`. Harmless today because every caller passes `t_`,
   but worth fixing while in the neighborhood.

## Risk Assessment

- **Phase 1**: very low. Additive, header-only, no existing call site changes behavior. The
  only failure mode is a compile error if `sieMin` is requested from a `UnitSystem` over a
  base EOS that lacks it.
- **Phase 2**: low but broad. Additive to the public API, but it adds two methods every
  model inherits, and the `ShiftedEOS` transform is a behavior *correction* relative to
  naive forwarding — it must land with the tests that pin it down.

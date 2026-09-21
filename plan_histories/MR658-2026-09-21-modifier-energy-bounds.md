# Plan: Expose Internal-Energy Bounds Through Modifiers (starting with `UnitSystem`)

**MR**: #XXX (rename this file once the MR number is assigned)
**Date**: 2026-09-19
**Status**: ✅ Phase 1 implemented — ✅ Phase 2 implemented (2026-09-21, at reviewer request)

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

#### Phase 2 Implementation Record (2026-09-21)

Carried out at a reviewer's request. Decisions taken on the open questions:

- **Naming**: `MinimumInternalEnergy()` / `MaximumInternalEnergy()`, matching the
  spelled-out style of the rest of the introspection API.
- **Concrete model scope**: all four models that can know their bounds —
  `SpinerEOSDependsRhoSie`, `StellarCollapse`, `EOSPAC`, and `SpinerEOSDependsRhoT`.
  `Helmholtz` keeps the permissive base defaults.
- **Negative `ScaledEOS` scale**: swap min/max when `scale_ < 0`, so the reported
  bounds stay ordered. The pre-existing `MinimumDensity`/`MaximumDensity` issue was
  deliberately *not* touched, to keep the diff scoped.

Deviation from the plan: step 2 proposed adding the two methods to
`SG_ADD_MODIFIER_INTROSPECTION_METHODS`, but `ShiftedEOS` uses that macro and needs a
real transform, so the macro would collide with its override. Instead the energy bounds
live in a second macro, `SG_ADD_MODIFIER_ENERGY_BOUNDS_METHODS`, applied only by the
modifiers that genuinely forward verbatim.

Files changed:

- `singularity-eos/eos/eos_base.hpp`: `MinimumInternalEnergy()`/`MaximumInternalEnergy()`
  defaults (`-1e100` / `1e100`); new `SG_ADD_MODIFIER_ENERGY_BOUNDS_METHODS` macro; fixed
  `SG_ADD_MODIFIER_INTROSPECTION_METHODS`'s parameter name (`t` → `t_`, open question 5).
- Verbatim forwarding via the new macro: `floored_energy.hpp`, `relativistic_eos.hpp`,
  `zsplit_eos.hpp` (each with a comment on why forwarding is right for that modifier).
- `ramps_eos.hpp`: hand-written forwarding, since it hand-writes its whole introspection
  block rather than using the macro.
- `shifted_eos.hpp`: explicit override adding `shift_` to both bounds.
- `scaled_eos.hpp`: explicit override multiplying by `scale_`, with the negative-scale swap.
- `eos_unitsystem.hpp`: explicit override dividing by the energy unit; the Phase 1
  `sieMin`/`sieMax` pass-throughs are kept, with a comment steering callers to the new
  methods.
- `eos_spiner_rho_sie.hpp`, `eos_stellar_collapse.hpp`: one-line delegation to the
  existing `sieMin()`/`sieMax()`.
- `eos_eospac.hpp`: new `sie_min_`/`sie_max_` members, filled from the `SesameMetadata`
  the constructor already fetches (answering open question 4: yes); ordering check added
  to `CheckParams`.
- `eos_spiner_rho_temp.hpp`: new `sie_min_`/`sie_max_` members plus a `setEnergyBounds_()`
  helper called from both construction paths (`loadDataboxes_` and the from-EOS
  constructor). Answering open question 3: the bounds are the extrema of the tabulated
  `sie_` field, unioned with `sieCold_`, *not* corner evaluations — corners would assume
  monotonicity in both arguments, which the table does not guarantee. Cached at load time
  because a `DataBox::min()` scan per call would be O(numRho*numT). `sieMin()`/`sieMax()`
  accessors added here too, for parity with the rho-sie table.
- `eos_variant.hpp`: two `visit`-based accessors.
- `python/module.hpp`: bound on the generic `eos_class<T>` template, so every bound type
  including the variant gets them.
- `test/test_eos_modifiers.cpp`: `BoundedGas` gains the two methods; new `THEN` blocks for
  the shifted+scaled composition, the negative-scale swap, relativistic pass-through, the
  unit-system conversion, and the permissive analytic defaults; new scenario exercising the
  motivating case — energy bounds read off a `singularity::Variant` holding a
  `UnitSystem<BoundedGas>`.
- `doc/sphinx/src/using-eos.rst`: documented both methods next to the density bounds, with
  a note distinguishing them from `MinInternalEnergyFromDensity` and an extended warning
  covering the negative default minimum.
- `doc/sphinx/src/modifiers.rst`: new "How Modifiers Transform the Energy Bounds" section
  enumerating each modifier's transform; the `UnitSystem` section now leads with the new
  API.
- `CHANGELOG.md`: second entry under `## Current develop` → `### Added`.
- Copyright years bumped and generative-AI notices added where missing.

Verification: `clang-format` (v21.1.4) applied to all changed C++ files. **The build and
test run were handed off to the user** — not verified by me.

#### Phase 1 Reversal: the `UnitSystem` pass-throughs were removed (2026-09-21)

A reviewer objected to the Phase 1 `sieMin()`/`sieMax()` pass-throughs in
`eos_unitsystem.hpp` — specifically to the lazy-instantiation construct, on the grounds
that the compile error a new developer hits is confusing even with the explanatory
comment. Their suggestion was `if constexpr` + a `static_assert`.

Resolution: **delete the pass-throughs instead.** Phase 2 makes them redundant —
`UnitSystem::MinimumInternalEnergy()` returns the same values, applies the same
`inv_sie_unit_` conversion, is unconditionally well formed for every `T` because `EosBase`
supplies a default, and works through the variant. The lazy-instantiation trick only
existed because `sieMin` was not part of the contract, which is exactly what Phase 2
fixed. This removes the construct rather than improving its diagnostic. Safe to do because
Phase 1 was never merged — `c4af3a91` lives only on `buechler/table_bounds`, and the only
in-tree caller of the pass-through was the Phase 1 test.

Also removed: the `sieMin()`/`sieMax()` aliases Phase 2 had added to
`SpinerEOSDependsRhoT` "for parity". Adding new instances of a name being signposted as
superseded works against the deprecation, so that table exposes only
`MinimumInternalEnergy`/`MaximumInternalEnergy`.

Deprecation of the *base-table* `sieMin`/`sieMax` (on `SpinerEOSDependsRhoSie` and
`StellarCollapse`) was scoped to **docs + CHANGELOG only** — a note in `models.rst` and an
entry under `### Deprecated`. No `[[deprecated]]` attribute. Two reasons:

1. **CI.** `SINGULARITY_STRICT_WARNINGS=ON` sets `-Wall -Werror` (`CMakeLists.txt:668`) and
   both `.github/workflows/warnings.yml` and `sanitizer.yml` enable it.
   `-Wdeprecated-declarations` is in `-Wall`, so the attribute turns every surviving
   in-tree call into a hard error: `python/module.cpp:143-144,157-158` (member pointers
   warn at bind time), `test/profile_stellar_collapse.cpp:121,146`, and the Phase 2
   delegations in `eos_spiner_rho_sie.hpp`/`eos_stellar_collapse.hpp`. There is also no
   C++ deprecation macro in the repo; the only precedent is the hand-rolled per-compiler
   `DEPRECATED_MODULE` in `singularity_eos.f90` (PR644), and these are
   `PORTABLE_FORCEINLINE_FUNCTION`, so an `ALLOW_DEPRECATED`-style escape hatch would be
   needed for device compilers.
2. **Family consistency.** `sieMin`/`sieMax` are one pair among eight table accessors.
   `rhoMin`/`rhoMax`/`TMin` already duplicate `MinimumDensity`/`MaximumDensity`/
   `MinimumTemperature` and have never been deprecated; `TMax`, `YeMin` and `YeMax` have
   **no** contract equivalent at all — there is no `MaximumTemperature()` in the library.
   `models.rst:2277` presents all eight as one set, so attributing only the energy pair
   would look arbitrary. Deprecating the family properly requires first adding
   `MaximumTemperature()` and deciding about `YeMin`/`YeMax` — a follow-on MR.

Python note: the new methods were bound on the generic `eos_class<T>` template, so Python
users already have the replacement. Giving them a real runtime `DeprecationWarning` would
mean replacing `def_property_readonly("sieMin", &T::sieMin)` with a lambda calling
`PyErr_WarnEx`; deferred with the rest of the attribute work.

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

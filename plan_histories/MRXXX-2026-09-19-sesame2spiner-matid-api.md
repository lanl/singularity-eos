# Add matid-based `saveAllMaterials` overloads to sesame2spiner-lib

## Context

A host code needs to hand sesame2spiner a list of sesame matids (idnos) and get an sp5
file back, without writing `.dat` parameter decks to disk first. Today the only entry
point is `saveAllMaterials(savename, filenames, ...)`
(`sesame2spiner/sesame2spiner/generate_files.cpp:173`), which parses config files as its
first step.

A colleague proposed adding two `vector<int>` overloads by copying the existing function
body. The signatures are right, but the copy duplicates ~90 lines that will drift from the
file-based path, and it promotes three latent defects into a public API. This plan keeps
the proposed signatures, makes the `(matids, params)` form the *single* implementation that
everything else delegates to, and splits out a core taking an already-open `hid_t` so the
host can also build one file up one material at a time.

Two decisions confirmed with the user:
- **Bad matid** → warn naming the matid, skip it, keep going, return nonzero. The good
  materials still land in the file.
- **Append** → expose the `hid_t` core plus root-attribute helpers. No new RAII class.

### Why a bad matid needs explicit handling

`eosGetMetadata` **discards** the `eosSafeLoad` return code
(`eospac-wrapper/eospac_wrapper.cpp:53`), and its `infoVals` buffers are zero-initialized
just above. A matid absent from the sesame file therefore returns metadata with
`rhoMin == rhoMax == 0` and no error of any kind, which flows into `getMatBounds` →
`log(0)` → NaN bounds. Tolerable when a human wrote the deck and read the spew;
a landmine when a host passes a programmatic list.

### Why the `log_type` root attribute needs validating

`SpinerEOS` hard-checks the **file-level** `log_type` attribute at load
(`singularity-eos/eos/eos_spiner_rho_temp.hpp:398-402`, `PORTABLE_ALWAYS_REQUIRE`). That
check cannot catch a file whose materials were written by builds with *different*
`log_type` settings: one material's grids will disagree with the attribute and be silently
mis-interpolated. Only matters once appending is possible, which is what this change
enables.

## Changes

### 1. `sesame2spiner/sesame2spiner/generate_files.hpp`

Replace the single `saveAllMaterials` declaration (line 59) with:

```cpp
// Root attributes (singularity_version, log_type). Call once on a newly created
// file before adding materials. Exposed for callers that own the hid_t themselves.
herr_t writeSP5RootAttributes(hid_t file);
// Verify an existing file's log_type matches this build. Called by the hid_t core.
herr_t checkSP5RootAttributes(hid_t file);

// Core. Adds materials to an already-open sp5 file. Safe to call repeatedly on the
// same file to build it up incrementally; duplicate matid and name detection read
// the file, so they stay correct across calls.
herr_t saveAllMaterials(hid_t file, const std::vector<int> &matids,
                        const std::vector<Params> &params, bool printMetadata,
                        Verbosity eospacWarn);

// Per-material parameter overrides, into a freshly created file.
herr_t saveAllMaterials(const std::string &savename, const std::vector<int> &matids,
                        const std::vector<Params> &params, bool printMetadata,
                        Verbosity eospacWarn);

// Standard sesame2spiner defaults for every material.
herr_t saveAllMaterials(const std::string &savename, const std::vector<int> &matids,
                        bool printMetadata, Verbosity eospacWarn);

// Existing file-based interface. Signature unchanged.
herr_t saveAllMaterials(const std::string &savename,
                        const std::vector<std::string> &filenames, bool printMetadata,
                        Verbosity eospacWarn);
```

Also: declare `bool checkMetadataValid(int matid, const SesameMetadata &metadata);` next to
the existing `checkValInMatBounds` (line 76), and drop the unused leading `int i` parameter
from `getMatBounds` (line 68).

### 2. `sesame2spiner/sesame2spiner/generate_files.cpp`

**The `hid_t` core is the only implementation.** It is the current loop body (lines
196-253) with these changes:

- Guard `matids.size() != params.size()` and `file < 0`; call `checkSP5RootAttributes`
  and bail on mismatch.
- **Duplicate detection via the file, not local state.** Drop `used_matids` /
  `used_names`; use `H5Lexists(file, ..., H5P_DEFAULT) > 0`. Both the `/<matid>` group and
  the `/<name>` soft link live at root (`saveMaterial` creates them at `loc`, lines 66-67),
  so one call covers each. Name collisions resolve by incrementing a suffix from `_2`
  until free — same numbering the current `used_names` map produces.
- **Validate metadata** with `checkMetadataValid` before use; on failure log the matid,
  increment a failure count, `continue`.
- **Per-material status.** Replace `status += saveMaterial(...)` plus the sticky
  `if (status != H5_SUCCESS)` (lines 248-252) with a local `herr_t` per material. Today,
  once any material fails, every *subsequent* material also prints "problem with HDF5"
  even on success, and the message never names a matid.
- **Count failures instead of summing `herr_t`.** `status += ...` lets a `-1` and a `+1`
  cancel to `0`, which the caller reads as success. Return `H5_SUCCESS` or `-1` based on a
  failure counter, preserving `main.cpp:59`'s `(status == H5_SUCCESS) ? 0 : 1`.
- **Pass `eospacWarn` to `eosGetMetadata`** instead of the hardcoded `Verbosity::Debug`
  (line 219). *This quiets the CLI's default output* — intended, but it is a visible
  change.
- Default name uses the matid, not the loop index: `"material_" + std::to_string(matid)`
  replaces `"material_" + std::to_string(i)` (line 209). The index is meaningless to a
  caller who passed matids.

**The three `savename` overloads become thin wrappers.** The file-based one keeps its
current parse loop (lines 176-182) and then delegates — `AddMaterials`
(`parser.hpp:91`) already produces exactly the `(params, matids)` pair the core wants:

```cpp
herr_t saveAllMaterials(const std::string &savename,
                        const std::vector<std::string> &filenames, bool printMetadata,
                        Verbosity eospacWarn) {
  std::vector<Params> params;
  std::vector<int> matids;
  for (auto const &filename : filenames) AddMaterials(params, matids, filename);
  return saveAllMaterials(savename, matids, params, printMetadata, eospacWarn);
}

herr_t saveAllMaterials(const std::string &savename, const std::vector<int> &matids,
                        bool printMetadata, Verbosity eospacWarn) {
  return saveAllMaterials(savename, matids, std::vector<Params>(matids.size()),
                          printMetadata, eospacWarn);
}
```

`Params` default-constructs (`parser.hpp:71`) and every `Get` call in the material loop
supplies a default, so an empty `Params` is exactly "standard defaults". The
`savename` + `(matids, params)` overload owns `H5Fcreate(H5F_ACC_TRUNC)` →
`writeSP5RootAttributes` → core → `H5Fclose`.

Finally, drop the `<unordered_map>` / `<unordered_set>` includes, now unused.

### 3. Docs

- `sesame2spiner/README.md` — currently a pure CLI transcript with no library section at
  all. Add "Using sesame2spiner as a library": the one-shot `vector<int>` call and the
  incremental `hid_t` loop, which is the host-facing documentation this whole change is
  for.
- `INTEGRATION_WITH_SESAME2SPINER_REFACTOR.md` — the de facto API doc. Update its
  `saveAllMaterials`/`getMatBounds` references (lines ~27, 34, 41, 125, 251, 268, 334-338).
- `CHANGELOG.md` — one entry.

`doc/sphinx/src/models.rst` needs no change: its sesame2spiner section (~line 1965) covers
the CLI and input-deck keys, neither of which changes.

## Deliberately out of scope

- **No `SP5File` RAII class.** Per the decision above; `writeSP5RootAttributes` plus the
  core's validation covers the hazard without a new public type.
- **`cmake/Format.cmake:56` has a typo** — `sesame2smpiner` instead of `sesame2spiner` —
  so `make format` silently skips all sesame2spiner `.cpp` files. I will *not* fix it here:
  it would reformat three files that have never been formatted and bury this diff. I'll run
  clang-format on just the edited regions and flag the typo separately.
- **Wiring up `sesame2spiner/test/test.cpp`.** It is registered by no CMake target
  anywhere, and is stale twice over: Catch2 v2 (`#define CATCH_CONFIG_MAIN` +
  `"catch.hpp"`) against the repo's Catch2 v3.7.1, and unqualified `eosGetMetadata` /
  `saveMaterial` / `eosDataOfRhoT` calls predating the `sesame2spiner` namespace. Porting
  it is real, unrelated work. **Caveat: this means the change lands with no automated
  coverage** — verification below is manual. `SINGULARITY_TEST_SESAME` (top-level
  `CMakeLists.txt:111-113`) is the correct gate if you want me to do it as a follow-up.

## Verification

Build with your usual kessel/spack flow; `SINGULARITY_BUILD_SESAME2SPINER` requires
`SINGULARITY_USE_SPINER`, `SINGULARITY_USE_SPINER_WITH_HDF5`, and
`SINGULARITY_USE_EOSPAC`.

1. **CLI regression — the critical one**, since the file-based path now delegates through
   new code. Capture `h5ls -r` output from the current binary *before* the change, then
   after: `./sesame2spiner -p examples/air.dat examples/steel.dat` must produce an
   identical file structure. Expect only stdout differences from the `eosGetMetadata`
   verbosity fix.
2. **Dedupe regression.** `./sesame2spiner examples/duplicate-test/*.dat` — `air.dat` and
   `air2.dat` are byte-identical, so this fixture exercises duplicate-matid skipping.
   Confirm the `H5Lexists` rewrite still skips it and yields the same structure as before.
3. **New one-shot API.** Short driver linking `sesame2spiner::sesame2spiner`:
   `saveAllMaterials("mats.sp5", {5030, 4272}, false, Verbosity::Quiet)`. `h5ls -r` must
   match what deck-driven run #1 produced.
4. **Bad matid.** `saveAllMaterials("mats.sp5", {5030, 99999, 4272}, ...)` — expect a
   stderr line naming 99999, a nonzero return, and an sp5 containing 5030 and 4272 with
   no NaN-bounded garbage group.
5. **Incremental path.** `H5Fcreate` + `writeSP5RootAttributes`, then
   `saveAllMaterials(file, {matid}, ...)` per matid, then `H5Fclose`. Structure must match
   the one-shot file from #3. Also call it twice with the same matid to confirm
   cross-call duplicate detection fires.
6. **Round-trip.** Load the incrementally built file through `SpinerEOSDependsRhoT` and
   evaluate a point or two, confirming it is a valid sp5 and the `log_type` check passes.

---

# Implementation Outcome (2026-09-19)

Implemented as planned, with one addition the plan did not anticipate.

## Files changed

- `sesame2spiner/sesame2spiner/generate_files.hpp` — new declarations
- `sesame2spiner/sesame2spiner/generate_files.cpp` — core + wrappers, fixes
- `sesame2spiner/README.md` — new "Using sesame2spiner as a library" section
- `CHANGELOG.md` — Added / Fixed / Changed entries
- `INTEGRATION_WITH_SESAME2SPINER_REFACTOR.md` — new "Programmatic saveAllMaterials
  API" section; annotated the older `SpinerTableGridParams` proposal as partially
  delivered; removed stale line-number references

No CMake changes were needed — no new translation units.

## Unplanned addition: overload ambiguity on braced string-literal lists

The plan did not catch this; a compile check of all call forms did.

```cpp
saveAllMaterials(savename, {"air.dat", "steel.dat"}, false, warn);  // ambiguous!
```

`std::vector<int>` has an iterator-pair constructor `vector(InputIt, InputIt)`, and
two `const char *` satisfy it: `iterator_traits<const char *>::value_type` is `char`,
which converts to `int`. So the braced list is viable as both a
`std::vector<std::string>` and a (nonsensical) `std::vector<int>`, and the call does
not compile.

Fixed by adding an inline overload taking `std::initializer_list<const char *>`,
which is an exact match for a braced list of string literals and therefore outranks
both `vector` candidates. It forwards to the `std::vector<std::string>` overload.

Note this only ever affected braced-literal call sites. `src/main.cpp` passes a named
`std::vector<std::string>`, so the CLI was never affected. The failure mode was a
compile error, not silent misbehavior.

## Behavior changes to be aware of

1. **CLI output is quieter.** `eosGetMetadata` now receives the caller's verbosity
   instead of a hardcoded `Verbosity::Debug`.
2. **Bad matids no longer produce garbage.** Previously an absent matid yielded
   all-zero metadata bounds, reaching `log()` and generating NaN grids with no
   diagnostic. Now reported by matid and skipped; other materials still save and the
   return is non-zero.
3. **Duplicate/name detection reads the file** (`H5Lexists`) rather than in-memory
   sets. Equivalent for a single batch call, and correct across incremental calls.
4. **`getMatBounds` lost its leading `int i`.** Nothing outside sesame2spiner called
   it.
5. **Default material names use the matid**, not the loop index: `material_<matid>`.

## Verification status

**Compiled, not run.** `build/`'s dependencies resolved into
`/tmp/buechler-ci-envs/spack`, which has since been cleaned, so `cmake --build`
cannot reconfigure. Compilation was verified by taking the exact flags from
`build/compile_commands.json` and substituting the in-tree `utils/ports-of-call` and
`utils/spiner` submodules for the missing spack prefixes. Clean under
`-Wall -Wextra` for `generate_files.cpp`, `src/main.cpp`, and a scratch driver
exercising all six call forms. Both edited files are clang-format clean (the
pre-existing code in them already was).

**Outstanding — needs the kessel/spack env and EOSPAC tables.** All six steps in the
Verification section above, most importantly:
- before/after `h5ls -r` comparison for
  `./sesame2spiner -p examples/air.dat examples/steel.dat`, since the file-based path
  now routes through new code
- `examples/duplicate-test/*.dat` to confirm the `H5Lexists` rewrite still skips the
  duplicate air matid
- a bad matid mixed into a real list
- incremental `hid_t` build compared against the equivalent one-shot file

## Known gaps

- **No automated coverage.** `sesame2spiner/test/test.cpp` is registered by no CMake
  target and would not compile if it were: Catch2 v2 idioms against the repo's v3.7.1,
  and unqualified calls predating the `sesame2spiner` namespace. Porting it under
  `SINGULARITY_TEST_SESAME` would give the bad-matid and cross-call dedupe paths real
  tests.
- **`cmake/Format.cmake:56` typo** — `sesame2smpiner` instead of `sesame2spiner` — so
  `make format` silently skips every sesame2spiner `.cpp`. Left alone deliberately;
  fixing it reformats files unrelated to this change. clang-format was run directly on
  the two edited files instead.

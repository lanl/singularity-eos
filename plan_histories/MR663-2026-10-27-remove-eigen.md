# Remove Eigen and Select Kokkos Kernels Automatically for Closures

**Date recorded:** 2026-10-27  
**Status:** Implemented; available CPU tests and package checks passed  
**Implementation commit:** `e9cf5eb6` (`remove eigen`)  
**MR:** Not assigned in this conversation

## Request and constraints

Remove the dependency on Eigen. Use Kokkos Kernels when closure solvers
require linear algebra, while keeping Kokkos and Kokkos Kernels optional
when mixed cell closures are disabled.

The user explicitly requested that we **not extend the regression suite**.
Validation was to use the existing regressions.

After the initial implementation, the user requested removal of the separate
Kokkos Kernels flag: enable the dependency automatically when closures need
it and Kokkos is enabled. The actual flag in the repository was named
`SINGULARITY_USE_KOKKOSKERNELS` (without the extra underscore in
`KOKKOS_KERNELS`).

## Findings that informed the plan

- The dense solve in `singularity-eos/closure/mixed_cell_models.hpp` already
  had a Kokkos Kernels implementation using KokkosBatched QR, ApplyQ, and Trsm.
  Eigen supplied an alternative CPU implementation.
- Eigen used column-pivoted Householder QR; the existing Kokkos Kernels path
  uses unpivoted QR. Switching CPU builds to the latter warranted running
  the existing numerical regressions.
- Three `kinetic_phasetransition_*` headers included linear algebra headers
  that they did not need.
- Closure headers and most closure tests were already conditionally built,
  but the `pte_2mat` example and Fortran test/install rules needed additional
  gating for closure-disabled builds.
- The original defaults enabled closures but disabled Kokkos. Removing the
  Eigen fallback required resolving that inconsistent combination.

## Original proposed plan

1. Enforce that closures require Kokkos and Kokkos Kernels in both submodule
   and external-package builds. Default closures to OFF, retaining Kokkos's
   OFF default. Initially retain the separate Kernels option for compatibility,
   but reject closure builds that disable it.
2. Remove the Eigen solver branch and unnecessary includes. Preserve the
   specialized 1x1 and 2x2 solves, solver interfaces, and scratch-space contract.
3. Remove Eigen discovery, linkage, CMake helpers, installed-package dependency
   declarations, submodule registration, and obsolete archive exclusions.
4. Audit closure-disabled builds, including examples, tests, headers, and the
   Fortran backend. Keep ordinary EOS functionality independent of the closure
   linear algebra dependency.
5. Validate numerical behavior, valid and invalid dependency combinations,
   submodule builds, and installed-package consumers. The original proposal
   included additional numerical coverage; the user explicitly removed that
   part before implementation. No regression cases were added or extended.
6. Update CI, presets, README/Sphinx documentation, and migration guidance.
   Note the corresponding changes needed in externally maintained Spack recipes.

## What was implemented

### Solver and headers

- Removed Eigen from `solve_Ax_b_wscr` and retained the existing KokkosBatched
  implementation for systems larger than 2x2.
- Preserved row scaling, specialized small-system solves, solver interfaces,
  scratch sizing, and the existing finite-result check.
- Guarded `mixed_cell_models.hpp` with `SINGULARITY_BUILD_CLOSURE` so disabled
  closures do not pull in the solver dependencies.
- Retained a compile-time diagnostic for attempting to build closures without
  the Kokkos portability strategy.
- Removed unused Eigen/Kokkos linear algebra includes from the three kinetic
  phase-transition headers.
- Updated the existing solver-test scratch-space comment; test behavior and
  regression coverage were not extended.

### Dependencies and configuration

- Changed `SINGULARITY_BUILD_CLOSURE` to default OFF.
- Added a configure-time error for closures ON with Kokkos OFF.
- Removed `SINGULARITY_USE_EIGEN`, its import/find/link logic,
  `cmake/singularity-eos/Eigen3.cmake`, and the Eigen branch in
  `cmake/submodule_configs.cmake`.
- Removed the `utils/eigen` submodule and its `.gitmodules` entry.
- Removed Eigen from the installed package's dependency discovery and from
  archive-script exclusions.
- Following the user's second request, removed the
  `SINGULARITY_USE_KOKKOSKERNELS` CMake option and compile definition entirely.
  Dependency import, discovery, linkage, and installed-package discovery now
  follow `SINGULARITY_BUILD_CLOSURE` directly.

The final contract is:

| Closures | Kokkos | Kokkos Kernels |
|---|---|---|
| OFF | OFF | Not required or selected by singularity-eos |
| OFF | ON | Not required or selected by singularity-eos |
| ON | OFF | Configuration error |
| ON | ON | Automatically required and linked |

Enable closures with:

```sh
cmake -S . -B build \
  -DSINGULARITY_BUILD_CLOSURE=ON \
  -DSINGULARITY_USE_KOKKOS=ON
```

Other options and dependency locations still depend on the chosen build mode.
No separate Kokkos Kernels enable flag is needed. External-package builds must
provide Kokkos Kernels; submodule builds import the bundled dependency.

### Build integration and documentation

- Gated `pte_2mat` on closures as well as its existing Spiner/HDF5 requirement.
- Gated the Fortran interface test on closures, and Fortran module installation
  on the existence of the compiled library target.
- Updated GitHub workflows to select closure configurations explicitly.
  The existing minimal Kokkos workflow now runs with closures ON and OFF;
  it uses the existing tests in both configurations.
- Updated Kokkos presets to enable closures explicitly.
- Removed references to the Kernels enable flag from presets, CI, and docs.
- Updated `README.build.md`, `README.spack.md`, and
  `doc/sphinx/src/building.rst` with the dependency contract and migration notes.
- Documented that external Spack recipes must require Kokkos and Kokkos Kernels
  for closures and stop depending on Eigen. Those external recipes were not
  modified in this workspace.

## Validation actually performed

The existing suites were rerun after removing the Kernels flag:

| Configuration | Existing tests | Result |
|---|---:|---|
| Closures ON, Kokkos Serial, automatic Kokkos Kernels; Spiner/HDF5, Helmholtz, and Fortran enabled | 107 | All passed |
| Closures OFF, Kokkos ON, no Kernels; Spiner enabled | 95 | All passed |
| Closures OFF, Kokkos OFF, no Kernels; Spiner disabled | 74 | All passed |

The enabled closure suite included the existing solver, PTE, and KPT unit
coverage and available PTE/Fortran regressions. No test failures were observed
in these configurations.

Additional checks:

- Built examples in the tested submodule configurations.
- Verified that closures ON with Kokkos OFF produces the intended error.
- Verified that an external-package closure build fails when Kokkos Kernels
  package discovery is disabled.
- Configured and installed external-package builds with closures both OFF
  and ON, with the removed Kernels option absent from their caches.
- Built and ran temporary installed-package consumers for an ordinary EOS
  and a closure solve. For the minimal consumer, Kokkos and Kokkos Kernels
  discovery were explicitly disabled. Eigen discovery was disabled for both.
- Verified the removed option was absent from the three regression-build
  caches and that no active repository references to the flag remained.
- Ran `git diff --check`, parsed `CMakePresets.json`, and ran the installed-header
  registration check successfully.

Temporary consumer checks were created under `/tmp`; they were not added to
the repository's regression suite. Build and test logs were also stored under
`/tmp`, including:

- `/tmp/singularity-noeigen-cpu-test.log`
- `/tmp/singularity-noeigen-kokkos-only-test.log`
- `/tmp/singularity-noeigen-minimal-test.log`

## Validation limitations and environment notes

- GPU builds/regressions were not run: CUDA and HIP tooling were unavailable.
- EOSPAC-dependent regressions were not run because that dependency was not
  available. The reported counts describe the enabled local configurations,
  not every optional regression in the repository.
- Local testing used GNU 15.2 and Kokkos/Kokkos Kernels 5.1.0 on aarch64.
- The first fresh configure needed network access to fetch Catch2. A subsequent
  configure encountered an upstream Kokkos Kernels cache-spelling issue;
  local configuration removed `KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV` and set
  `KokkosKernels_ENABLE_SUPERNODAL_SPTRSV=OFF`. No upstream source fix was made.
- Existing unrelated working-tree edits were preserved during implementation.
  At the time this history was written, the removal changes were committed as
  `e9cf5eb6`; the pre-existing untracked content in `utils/kokkos-kernels`
  remained untouched.

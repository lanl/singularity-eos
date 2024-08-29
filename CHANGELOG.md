# Changelog

## Current develop

### Added (new features/APIs/variables/...)
- [[PR361]](https://github.com/lanl/singularity-eos/pull/361) Added tests for PTEsolver and added a strawman kinetic phase transition framework
- [[PR330]](https://github.com/lanl/singularity-eos/pull/330) Piecewise grids for Spiner EOS.

### Fixed (Repair bugs, etc)
- [[PR330]](https://github.com/lanl/singularity-eos/pull/330) Includes a fix for extrapolation of specific internal energy in SpinerEOS.
- [[PR401]](https://github.com/lanl/singularity-eos/pull/401) Fix for internal energy scaling in PTE closure
- [[PR403]](https://github.com/lanl/singularity-eos/pull/403) Fix Gruneisen EOS DensityEnergyFromPressureTemperature function

### Changed (changing behavior/API/variables/...)
- [[PR407]](https://github.com/lanl/singularity-eos/pull/407) Update C++ standard to C++17

### Infrastructure (changes irrelevant to downstream codes)
- [[PR402]](https://github.com/lanl/singularity-eos/pull/402) Added stiff gas to python interface

### Removed (removing behavior/API/varaibles/...)

## Release 1.9.0
Date: 7/29/2024

### Added (new features/APIs/variables/...)
- [[PR377]](https://github.com/lanl/singularity-eos/pull/377) Moved much of the variant creating machinery and initialization machinery into separate header files. This is useful for downstream codes that use custom variants and helps with producing plugins.
- [[PR292]](https://github.com/lanl/singularity-eos/pull/292) Added Carnahan-Starling EoS
- [[PR#362]](https://github.com/lanl/singularity-eos/pull/362) Add lambda to thermalqs
- [[PR#339]](https://github.com/lanl/singularity-eos/pull/339) Added COMPONENTS to singularity-eos CMake install, allowing to select a minimal subset needed e.g. for Fortran bindings only
- [[PR#336]](https://github.com/lanl/singularity-eos/pull/336) Included code and documentation for a full, temperature consistent, Mie-Gruneisen EOS based on a pressure power law expansion in eta = 1-V/V0. PowerMG.
- [[PR334]](https://github.com/lanl/singularity-eos/pull/334) Include plugins infrastructure
- [[PR331]](https://github.com/lanl/singularity-eos/pull/331) Included code and documentation for a full, temperature consistent, Mie-Gruneisen EOS based on a linear Us-up relation. MGUsup.
- [[PR326]](https://github.com/lanl/singularity-eos/pull/326) Document how to do a release
- [[PR#357]](https://github.com/lanl/singularity-eos/pull/357) Added support for C++17 (e.g., needed when using newer Kokkos).
- [[PR#382]](https://github.com/lanl/singularity-eos/pull/382) Added debug checks to the `get_sg_eos()` interface to ensure sane values are returned

### Fixed (Repair bugs, etc)
- [[PR380]](https://github.com/lanl/singularity-eos/pull/380) Set material internal energy to 0 if not participating in the pte solve to make sure potentially uninitialized data is set.
- [[PR370]](https://github.com/lanl/singularity-eos/pull/370) Fix bulk modulus calculation in spiner EOS
- [[PR343]](https://github.com/lanl/singularity-eos/pull/343) Add chemical potentials to stellar collapse gold files
- [[PR342]](https://github.com/lanl/singularity-eos/pull/342) Fix missing using statement in stellar collapse root finding routines
- [[PR341]](https://github.com/lanl/singularity-eos/pull/341) Short-circuit HDF5 machinery when cray-wrappers used in-tree
- [[PR340]](https://github.com/lanl/singularity-eos/pull/335) Fix in-tree builds with plugin infrastructure
- [[PR335]](https://github.com/lanl/singularity-eos/pull/335) Fix missing hermite.hpp in CMake install required for Helmholtz EOS
- [[PR356]](https://github.com/lanl/singularity-eos/pull/356) Guard against FPEs in the PTE solver
- [[PR356]](https://github.com/lanl/singularity-eos/pull/356) Update CMake for proper Kokkos linking in Fortran interface
- [[PR373]](https://github.com/lanl/singularity-eos/pull/373) Initialize cache in `get_sg_eos*` functions
- [[PR374]](https://github.com/lanl/singularity-eos/pull/374) Make the Davis EOS more numerically robust
- [[PR383]](https://github.com/lanl/singularity-eos/pull/383) Fix bug in step scaling for PTE solver

### Changed (changing behavior/API/variables/...)
- [[PR363]](https://github.com/lanl/singularity-eos/pull/363) Template lambda values for scalar calls
- [[PR372]](https://github.com/lanl/singularity-eos/pull/372) Removed E0 from Davis Products EOS in favor of using the shifted EOS modifier. CHANGES API!
- [[PR#382]](https://github.com/lanl/singularity-eos/pull/382) Changed `get_sg_eos()` API to allow optionally specifying the mass fraction cutoff for materials to participate in the PTE solver

### Infrastructure (changes irrelevant to downstream codes)
- [[PR329]](https://github.com/lanl/singularity-eos/pull/329) Move vinet tests into analytic test suite
- [[PR328]](https://github.com/lanl/singularity-eos/pull/328) Move to catch2 v3

### Removed (removing behavior/API/varaibles/...)

## Release 1.8.0
Date: 11/28/2023

### Fixed (Repair bugs, etc)
- [[PR278]](https://github.com/lanl/singularity-eos/pull/278) Fixed EOSPAC unit conversion errors for scalar lookups
- [[PR316]](https://github.com/lanl/singularity-eos/pull/316) removed `fmax-errors=3` from `singularity-eos` compile flags
- [[PR296]](https://github.com/lanl/singularity-eos/pull/296) changed `CMAKE_SOURCE_DIR` to `PROJECT_SOURCE_DIR` to fix downstream submodule build
- [[PR291]](https://github.com/lanl/singularity-eos/pull/291) package.py updates to reflect new CMake options
- [[PR290]](https://github.com/lanl/singularity-eos/pull/290) Added target guards on export config
- [[PR288]](https://github.com/lanl/singularity-eos/pull/288) Don't build tests that depend on spiner when spiner is disabled
- [[PR287]](https://github.com/lanl/singularity-eos/pull/287) Fix testing logic with new HDF5 options
- [[PR282]](https://github.com/lanl/singularity-eos/pull/282) Fix missing deep copy in sap polynomial tests
- [[PR281]](https://github.com/lanl/singularity-eos/pull/281) Pin spiner in spackage to a specific, tested version
- [[PR267]](https://github.com/lanl/singularity-eos/pull/267) Add missing eosSafeDestroy call in EOSPAC::Finalize
- [[PR263]](https://github.com/lanl/singularity-eos/pull/263) Fix stellar collapse mass fraction interpolation
- [[PR258]](https://github.com/lanl/singularity-eos/pull/258) Fix EOSPAC vector performance and add warnings when slower version gets selected.
- [[PR243]](https://github.com/lanl/singularity-eos/pull/243) Remove undefined behavior caused by diagnostic variables. Also fixed some compiler warnings.
- [[PR239]](https://github.com/lanl/singularity-eos/pull/239#issuecomment-1473166925) Add JQuery to dependencies for docs to repair build
- [[PR238]](https://github.com/lanl/singularity-eos/pull/238) Fixed broken examples
- [[PR228]](https://github.com/lanl/singularity-eos/pull/228) added untracked header files in cmake
- [[PR215]](https://github.com/lanl/singularity-eos/pull/215) and [[PR216]](https://github.com/lanl/singularity-eos/pull/216) fix duplicate definition of EPS and fix CI
- [[PR232]](https://github.com/lanl/singularity-eos/pull/228) Fixed uninitialized cmake path variables
- [[PR308]](https://github.com/lanl/singularity-eos/pull/308) spack builds +fortran now compile via correct blocking out of interfaces via preprocessor ifdef

### Added (new features/APIs/variables/...)
- [[PR338]](https://github.com/lanl/singularity-eos/pull/338) added chemical potentials from EoS
- [[PR269]](https://github.com/lanl/singularity-eos/pull/269) Add SAP Polynomial EoS
- [[PR278]](https://github.com/lanl/singularity-eos/pull/278) Added EOSPAC option functionality in class constructor
- [[PR278]](https://github.com/lanl/singularity-eos/pull/278) Added a new function for returning the minimum energy as a function of density for an EOS (only EOSPAC at the moment)
- [[PR278]](https://github.com/lanl/singularity-eos/pull/278) Added a new Fortran API for simple pressure and bulk moduli lookups
- [[PR306]](https://github.com/lanl/singularity-eos/pull/306) Added generic Evaluate method
- [[PR304]](https://github.com/lanl/singularity-eos/pull/304) added a Newton-Raphson root find for use with the Helmholtz EoS
- [[PR265]](https://github.com/lanl/singularity-eos/pull/265) Add missing UnitSystem modifier combinations to variant and EOSPAC
- [[PR279]](https://github.com/lanl/singularity-eos/pull/279) added noble-abel EoS
- [[PR274]](https://github.com/lanl/singularity-eos/pull/274) added a stiffened gas EoS
- [[PR264]](https://github.com/lanl/singularity-eos/pull/274) added a Helmholtz EoS
- [[PR254]](https://github.com/lanl/singularity-eos/pull/254) added ability to peel off modifiers as needed
- [[PR250]](https://github.com/lanl/singularity-eos/pull/250) added mass fraction output to stellar collapse eos
- [[PR226]](https://github.com/lanl/singularity-eos/pull/226) added entropy interpolation to stellar collapse eos
- [[PR202]](https://github.com/lanl/singularity-eos/pull/202) added the Vinet analytical EOS wth test cases and documentation.
- [[PR226]](https://github.com/lanl/singularity-eos/pull/226) added entropy interpolation to stellar collapse eos
- [[PR209]](https://github.com/lanl/singularity-eos/pull/209) added more documentation around how to contribute to the project and also what a contributor can expect from the core team
- [[PR214]](https://github.com/lanl/singularity-eos/pull/214) added documentation about adding a new EOS
- [[PR228]](https://github.com/lanl/singularity-eos/pull/228) and [[PR229]](https://github.com/lanl/singularity-eos/pull/229) added untracked header files in cmake
- [[PR233]](https://github.com/lanl/singularity-eos/pull/233) Added entropy infrastructure to singularity. Add entropy for the ideal gas and modifiers while throwing an error for EOS where entropy is not implemented yet
- [[PR177]](https://github.com/lanl/singularity-eos/pull/177) added EOSPAC vector functions

### Changed (changing behavior/API/variables/...)
- [[PR311]](https://github.com/lanl/singularity-eos/pull/311) Refactor EOSBuilder for flexibility and extensibility
- [[PR310]](https://github.com/lanl/singularity-eos/pull/310) Speed up and clean up tests
- [[PR295]](https://github.com/lanl/singularity-eos/pull/295) Add fast logs to singularity-eos
- [[PR246]](https://github.com/lanl/singularity-eos/pull/246) Update CMake with upstream changes
- [[PR223]](https://github.com/lanl/singularity-eos/pull/223) Update ports-of-call and add portable error handling
- [[PR234]](https://github.com/lanl/singularity-eos/pull/234) update ports-of-call to correct for undefined behavior in error handling
- [[PR219]](https://github.com/lanl/singularity-eos/pull/219) Removed static analysis from re-git pipeline
- [[PR233]](https://github.com/lanl/singularity-eos/pull/233) Exposed entropy for the EOS type (now required for future EOS)
- [[PR308]](https://github.com/lanl/singularity-eos/pull/308) Fortran initialization interface functions no longer require modifier arrays, they are optional parameters.

### Infrastructure (changes irrelevant to downstream codes)
- [[PR190]](https://github.com/lanl/singularity-eos/pull/190) update CI on re-git
- [[PR245]](https://github.com/lanl/singularity-eos/pull/245) Separating get_sg_eos to other files. Build/compilation improvements, warning fixes/suppression.
- [[PR308]](https://github.com/lanl/singularity-eos/pull/308) Added a fortran test.

### Removed (removing behavior/API/varaibles/...)
- [[PR293]](https://github.com/lanl/singularity-eos/pull/293) Removing PTofRE function. This will no longer be callable downstream.

## Release 1.7.0
Date: 12/14/2022

### Fixed (Repair bugs, etc)

### Added (new features/APIs/variables/...)
- [[PR209]](https://github.com/lanl/singularity-eos/pull/209) added more documentation around how to contribute to the project and also what a contributor can expect from the core team

### Changed (changing behavior/API/variables/...)
- [[PR192]](https://github.com/lanl/singularity-eos/pull/192) remove h5py and start using gold files

### Fixed (not changing behavior/API/variables/...)
- [[PR198]](https://github.com/lanl/singularity-eos/pull/198) Use a more robust HDF5 handling that eliminates some edge breaks
- [[PR181]](https://github.com/lanl/singularity-eos/pull/181) change from manual dependecy handling to using hdf5 interface targets

### Infrastructure (changes irrelevant to downstream codes)
- [[PR206]](https://github.com/lanl/singularity-eos/pull/206) add another way to build sphinx docs

### Removed (removing behavior/API/varaibles/...)

## Release 1.6.2
Date: 10/12/2022

### Fixed (Repair bugs, etc)
- [[PR183]](https://github.com/lanl/singularity-eos/pull/183) fortify cmake export config to always have interface targets of dependencies that need them
- [[PR174]](https://github.com/lanl/singularity-eos/pull/174) fix build configuration when closures are disabled
- [[PR157]](https://github.com/lanl/singularity-eos/pull/157) fix root finders in Gruneisen, Davis, and JWL
- [[PR151]](https://github.com/lanl/singularity-eos/pull/151) fix module install

### Added (new features/APIs/variables/...)
- [[PR175]](https://github.com/lanl/singularity-eos/pull/175) document some builds
- [[PR164]](https://github.com/lanl/singularity-eos/pull/164) provide facilities for an initial temperature guess for PTE
- [[PR156]](https://github.com/lanl/singularity-eos/pull/156) This PR adds 2 new PTE solvers, solvers that obtain a PTE solution when either P or T are known. This enables mixed material rho-T and rho-P initializations to be dealt with correctly.

### Changed (changing behavior/API/variables/...)
- [[PR156]](https://github.com/lanl/singularity-eos/pull/156) This PR changes how the get_sg_eos function calls PTE solvers and does lookups based on input condition. Each input condition now calls into its own specialized kernel that ensures the inputs and outputs are treated appropriately.
- [[PR168]](https://github.com/lanl/singularity-eos/pull/168) move EOS files to header-only
- [[PR167]](https://github.com/lanl/singularity-eos/pull/167) allow for the possiblity Kokkos version can't be inferred

## Release 1.6.1
Date: 07/07/2022

### Added (new features/APIs/variables/...)
- [[PR150]](https://github.com/lanl/singularity-eos/pull/150) This changelog

### Changed (changing behavior/API/variables/...)
- [[PR146]](https://github.com/lanl/singularity-eos/pull/146) Changes needed to cmake to enable spackage upstream

### Fixed (not changing behavior/API/variables/...)

### Infrastructure (changes irrelevant to downstream codes)

### Removed (removing behavior/API/varaibles/...)

## Release 1.6.0
Date: 07/07/2022

This is the start of changelog

© 2021-2024. Triad National Security, LLC. All rights reserved.  This
program was produced under U.S. Government contract 89233218CNA000001
for Los Alamos National Laboratory (LANL), which is operated by Triad
National Security, LLC for the U.S.  Department of Energy/National
Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of
Energy/National Nuclear Security Administration. The Government is
granted for itself and others acting on its behalf a nonexclusive,
paid-up, irrevocable worldwide license in this material to reproduce,
prepare derivative works, distribute copies to the public, perform
publicly and display publicly, and to permit others to do so.

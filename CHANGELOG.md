# Changelog

## Current develop
- [[PR269]](https://github.com/lanl/singularity-eos/pull/269) Add SAP Polynomial EoS

### Fixed (Repair bugs, etc)
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

### Added (new features/APIs/variables/...)
- [[PR274]](https://github.com/lanl/singularity-eos/pull/274) added a stiffened gas EoS
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
- [[PR246]](https://github.com/lanl/singularity-eos/pull/246) Update CMake with upstream changes
- [[PR223]](https://github.com/lanl/singularity-eos/pull/223) Update ports-of-call and add portable error handling
- [[PR234]](https://github.com/lanl/singularity-eos/pull/234) update ports-of-call to correct for undefined behavior in error handling
- [[PR219]](https://github.com/lanl/singularity-eos/pull/219) Removed static analysis from re-git pipeline
- [[PR233]](https://github.com/lanl/singularity-eos/pull/233) Exposed entropy for the EOS type (now required for future EOS)

### Infrastructure (changes irrelevant to downstream codes)
- [[PR190]](https://github.com/lanl/singularity-eos/pull/190) update CI on re-git

### Removed (removing behavior/API/varaibles/...)

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

© 2021-2023. Triad National Security, LLC. All rights reserved.  This
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

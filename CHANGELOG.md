# Changelog

## Current develop

### Fixed (Repair bugs, etc)
- [[PR174]](https://github.com/lanl/singularity-eos/pull/174) fix build configuration when closures are disabled
- [[PR157]](https://github.com/lanl/singularity-eos/pull/157) fix root finders in Gruneisen, Davis, and JWL
- [[PR151]](https://github.com/lanl/singularity-eos/pull/151) fix module install

### Added (new features/APIs/variables/...)
- [[PR182]](https://github.com/lanl/singularity-eos/pull/182) add fortran test
- [[PR175]](https://github.com/lanl/singularity-eos/pull/175) document some builds
- [[PR164]](https://github.com/lanl/singularity-eos/pull/164) provide facilities for an initial temperature guess for PTE
- [[PR156]](https://github.com/lanl/singularity-eos/pull/156) This PR adds 2 new PTE solvers, solvers that obtain a PTE solution when either P or T are known. This enables mixed material rho-T and rho-P initializations to be dealt with correctly.

### Changed (changing behavior/API/variables/...)
- [[PR156]](https://github.com/lanl/singularity-eos/pull/156) This PR changes how the get_sg_eos function calls PTE solvers and does lookups based on input condition. Each input condition now calls into its own specialized kernel that ensures the inputs and outputs are treated appropriately.
- [[PR168]](https://github.com/lanl/singularity-eos/pull/168) move EOS files to header-only
- [[PR167]](https://github.com/lanl/singularity-eos/pull/167) allow for the possiblity Kokkos version can't be inferred

### Fixed (not changing behavior/API/variables/...)

### Infrastructure (changes irrelevant to downstream codes)

### Removed (removing behavior/API/varaibles/...)

### Added (modifiers with python bindings)

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

© 2021-2022. Triad National Security, LLC. All rights reserved.  This
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

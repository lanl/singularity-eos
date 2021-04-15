Singularity EOS
===

[![Tests](https://github.com/lanl/singularity-eos/actions/workflows/tests.yml/badge.svg)](https://github.com/lanl/singularity-eos/actions/workflows/tests.yml)

A collection of closure models and tools useful for multiphysics codes

## Components

- `eos` contains the analytic and tabular equation of state infrastructure
- `mix` contains mixed cell closure models. Currently pressure-temperature equilibrium is supported. Requires `eos`.
- `utils/sesame2spiner` converts sesame tables to 

## To Build

At its most basic:
```bash
git clone --recursive git@gitlab.lanl.gov:singularity/singularity-eos.git
cd singularity-eos
mkdir bin
cd bin
cmake ..
make
```

### To Make tests

To build the tests, in the `singularity-eos` root directory, clone `singularity-eos-data`, which is a [git-lfs](https://git-lfs.github.com/) directory. First you need to install git-lfs. This only needs to be done once. Download git-lfs via, e.g.,
```bash
apt install git-lfs
```
and then add the requisite files into your home directory via
```bash
git lfs install
```

Then, in the singularity-eos root directory:
```bash
git clone git@gitlab.lanl.gov:singularity/singularity-eos-data.git
mkdir -p bin
cd bin
cmake -DSINGULARITY_BUILD_TESTS=ON ..
make -j
make test
```

### Build Options

A number of options are avaialable for compiling:

| Option                            | Default | Comment                                                                              |
| --------------------------------- | ------- | ------------------------------------------------------------------------------------ |
| SINGULARITY_USE_HDF5              | ON      | Enables HDF5. Required for SpinerEOS, StellarCollapseEOS, and sesame2spiner          |
| SINGULARITY_USE_FORTRAN           | ON      | Enable Fortran API for equation of state                                             |
| SINGULARITY_USE_KOKKOS            | OFF     | Uses Kokkos as the portability backend. Currently only Kokkos is supported for GPUs. |
| SINGULARITY_USE_EOSPAC            | OFF     | Link against EOSPAC. Needed for sesame2spiner and some tests.                        |
| SINGULARITY_USE_CUDA              | OFF     | Target nvidia GPUs via cuda. Currently requires Kokkos.                              |
| SINGULARITY_USE_KOKKOSKERNELS     | OFF     | Use Kokkos Kernels for linear algebra. Needed for mixed cell closure models on GPU   |
| SINGULARITY_BUILD_MIX             | ON      | Builds mixed cell closure machinery for multi-material problems                      |
| SINGULARITY_BUILD_TESTS           | OFF     | Build test infrastructure.                                                           |
| SINGULARITY_TEST_SESAME           | OFF     | Test the Sesame table readers                                                        |
| SINGULARITY_TEST_STELLAR_COLLAPSE | OFF     | Test the Stellar Collapse table readers                                              |
| SINGULARITY_BUILD_SESAME2SPINER   | OFF     | Builds the conversion tool sesame2spiner which makes files readable by SpinerEOS     |
| SINGULARITY_INVERT_AT_SETUP       | OFF     | For tests, pre-invert eospac tables.                                                 |
| SINGULARITY_BETTER_DEBUG_FLAGS    | ON      | Enables nicer GPU debug flags. May interfere with in-tree builds as a submodule      |
| SINGULARITY_HIDE_MORE_WARNINGS    | OFF     | Makes warnings less verbose. May interfere with in-tree builds as a submodule        |
| SINGULARITY_SUBMODULE_MODE        | OFF     | Sets cmake flags for in-tree builds                                                  |

### Building on Darwin power9 nodes with the IBM xl compiler

```bash
git clone --recursive https://gitlab.lanl.gov/singularity/singularity-eos.git
cd singularity-eos
mkdir bin && cd bin
source /usr/projects/eap/dotfiles/.bashrc xl reset
module switch ibm/xlc-16.1.1.5-xlf-16.1.1.5 ibm/xlc-16.1.1.7-xlf-16.1.1.7-gcc-7.4.0
export FC="xlf_r -F/usr/projects/opt/ppc64le/ibm/xlf-16.1.1.7/xlf/16.1.1/etc/xlf.cfg.rhel.7.8.gcc.7.4.0.cuda.10.1"
module load cmake/3.17.3
cd /build/directory
cmake -DCMAKE_BUILD_TYPE=Release -DSINGULARITY_USE_KOKKOS=ON -DSINGULARITY_USE_HDF5=ON -DSINGULARITY_USE_EOSPAC=ON -DSINGULARITY_EOSPAC_INSTALL_DIR=/usr/projects/eap/spack/spack-develop_20191010/opt/spack/linux-rhel7-power9le/xl_r-16.1.1.7-plain/eospac-6.4.0beta.2-mp6hd267l4iokwr6ght2kp5ukk5zc6vq -DSINGULARITY_BUILD_TESTS=ON ..
make -j
make test
```

## Copyright

© 2021. Triad National Security, LLC. All rights reserved.  This
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

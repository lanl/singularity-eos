These instructions document the steps to using Spack to set up a development environment for `singularity-eos`.
This documentation will assume as little familiarity with Spack as is possible.
However, there may be instances where familiarty with Spack and CMake are required to successfully deploy at the developer end.

## Preliminaries
### Clone the Spack repository
Spack is a package management project similar in use to software like `conan` or `conda`, with a focus for HPC. You can find the documentation [here](https://spack.readthedocs.io/en/latest/index.html), and the repository [here](https://github.com/spack/spack).

I clone it to `${HOME}/spack` for convience, but the location does not matter. Moving forward, I will refer to this directory as `${SPACK_TOP}` (this env var is only for reference, and does not gets set in the process)

```bash
$> cd
$> git clone https://github.com/spack/spack.git 
```

### Do initial spack configuration
Here we will generate a base spack configuration for your local machine. There are ways and there are ways to do this, but we will try to keep it as simple as possible.

##### Step 0: (optional)
In the following steps, Spack will try to look for compilers and packages it can automatically generate configurations for.
It can be somewhat clever in where it looks, but if you want it to find some particular package or compiler, you're best bet is to bring it into your `${PATH}`, for instance by doing a `module load X` command, e.g.

```bash
$> module load gcc/9.3.0 openmpi/4.0.5-gcc9.3.0 cmake/3.19.2 cuda/11.4
```
##### Step 1: activate spack
This step will bring in all the spack machinary into your path.
```bash
$> . ${SPACK_TOP}/share/spack/setup-env.sh
```
##### Step 2: Detect compilers and packages
This step will add a local configuration to `${HOME}/.spack`, specific to your workspace machine.
```bash
$> spack compiler find
$> spack external find --not-buildable
```
Each command will generate output telling you what was found, and where the configuration was placed.
Note that this process is not perfect. Some packages may be missed, or their specific versions will not be usable as a dependency for later concretization.
It is possible for spack to find no compilers or pre-installed packages, spack will just build everything (even compilers). This step, however, hopefully alievates some unnecessary work.

## Installing the Dependency package
This will guide you through using the `singularity-eos` dependency package to install the required dependencies.

##### Step 1: Add the `singularity-eos` Spack repository
Go to the top directory of the source code and add the repo,
```bash
$> cd <path/to/source/top>
$> spack repo add spack-repo/
```
##### Step 2: Use the `singularity-eos` Spack dependency package to install
The following variants may be used to customize the dependencies and their build variants

```bash
    Name [Default]          Allowed values          Description
    ====================    ====================    ========================================

    build_extra [none]      none,                   Build converters
                            stellarcollapse,        
                            sesame                  
    cuda [off]              on, off                 Build with CUDA
    cuda_arch [none]        none, 53, 35, 10,       CUDA architecture
                            61, 70, 20, 50, 21,     
                            37, 52, 72, 75, 12,     
                            60, 32, 11, 80, 13,     
                            86, 62, 30              
    enable_fortran [on]     on, off                 Enable building fortran interface
    enable_tests [off]      on, off                 Build tests
    kokkos [on]             on, off                 Enable kokkos
    kokkos-kernels [off]    on, off                 Enable kokkos-kernals for linear algebra
    mpi [off]               on, off                 Build with MPI support
    openmp [off]            on, off                 Enable openmp
```
These variants are subject to change as development continues. You may see all up-to-date variants using spack with
```bash
$> spack info singularity-eos-deps
```
The default installation, which installs very minimal dependencies
```bash
$> spack install singularity-eos-deps
```
##### Step 3: Load the package

Once installed, loading the package will load all the dependencies into your `${PATH}`
```bash
$>spack load singularity-eos-deps
```
## Building `singularity-eos` with Spack dependencies
Once loaded, you may run CMake as you like.

A useful feature of using the Spack dependencies package is automatic CMake configuration. Loading `singularity-eos-deps` will define an shell environment variable `${SINGULARITYEOS_TCF}`, which points to a CMake source with pre-configured CMake cache variables which were generated based on the selected variants.

For example
```bash
$> cd <build-dir>
$> cd <path/to/build/dir>
$> cmake -C $SINGULARITYEOS_TCF -DCMAKE_BUILD_TYPE=Debug <path/to/source>
$> cmake --build . -j
$> ctest
```
These pre-configured options may be overriden explicitly, e.g.
```bash
$> cmake -C $SINGULARITYEOS_TCF -DSINGULARITY_USE_FORTRAN=ON -DCMAKE_BUILD_TYPE=Debug <path/to/source>
```
The trailing options will override any set in `${SINGULARITYEOS_TCF}`.
## Cleanup

To unload, simply

`spack unload singularity-eos-deps`

which will return your shell to a state prior to loading.

## Some example variants
This non-exhaustive list provides examples for building with different back-end variants

To build with `kokkos` using a `cuda` backend (and required CUDA architecture)
```bash
$> spack install singularity-eos-deps+cuda+kokkos+kokkos-kernels cuda_arch=70 +enable_tests+enable_fortran build_extra=sesame,stellarcollapse
```

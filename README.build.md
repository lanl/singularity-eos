# Overview 

The `singularity-eos` build system is designed with two goals in mind

1. Portability to a wide range of host codes, system layouts, and underlying hardware
2. Ease of code development, and flexibility for developers

These considerations continue to guide development of the tools and workflows 
in working with `singularity-eos`.

## Basics

The build of `singularity-eos` can take two forms:

1. Submodule mode
2. Standalone mode

These will be described in more detail below, but in brief _submodule mode_ is intended 
for downstream codes that build `singularity-eos` source code directly in the build 
(sometimes referred to as "in-tree"), while _standalone mode_ will build `singularity-eos` 
as an independent library that can be installed onto the system.

The most important distinction between the modes is how dependencies are handled.
_submodule mode_ will use *internal* source clones of key dependencies (located in 
`utils\`), effectively building these dependencies as part of the overall `singularity-eos`
build procedure. It should be noted, however, that there are optional dependencies 
that are not provided internally and must be separately available.

In _standalone mode_, *all* dependencies must be available in the environment,
and be discoverable to CMake. While not required, it is encouraged to use the 
dependency management tool `spack` to help facilitate constructing a build environment, 
as well as deploying `singularity-eos`. Example uses of `spack` for these purposes 
are provided below.

A CMake configuration option is provided that allows developers to select a specific 
mode (`SINGULARITY_FORCE_SUBMODULE_MODE`), however this is intended for internal development 
only. The intended workflow is to let `singularity-eos` decide that appropriate mode, which 
it decides based on inspecting the project directory that the source resides in.

## Options for configuring the build

Most configuration options are the same between the two builds. _standalone_ / _submodule_ specific options are touched on in the sections detailing those build modes.

The main CMake options to configure building are in the following table: 

| Option | Default | Comment |
|--|--|--|
| `SINGULARITY_USE_SPINER`            | ON      | Enables EOS objects that use `spiner`.|
| `SINGULARITY_USE_FORTRAN`           | ON      | Enable Fortran API for equation of state.|
| `SINGULARITY_USE_KOKKOS`            | OFF     | Uses Kokkos as the portability backend. Currently only Kokkos is supported for GPUs.|
| `SINGULARITY_USE_EOSPAC`            | OFF     | Link against EOSPAC. Needed for sesame2spiner and some tests.|
| `SINGULARITY_BUILD_TESTS`           | OFF     | Build test infrastructure.|
| `SINGULARITY_BUILD_PYTHON`          | OFF     | Build Python bindings.|
| `SINGULARITY_INVERT_AT_SETUP`       | OFF     | For tests, pre-invert eospac tables.|
| `SINGULARITY_BETTER_DEBUG_FLAGS`    | ON      | Enables nicer GPU debug flags. May interfere with in-tree builds as a submodule.|
| `SINGULARITY_HIDE_MORE_WARNINGS`    | OFF     | Makes warnings less verbose. May interfere with in-tree builds as a submodule.|
| `SINGULARITY_FORCE_SUBMODULE_MODE`  | OFF     | Force build in _submodule_ mode.|
| `SINGULARITY_USE_SINGLE_LOGS`  | OFF     | Use single precision logarithms (may degrade accuracy).|
| `SINGULARITY_USE_TRUE_LOG_GRIDDING`  | OFF     | Use grids that conform to logarithmic spacing.|

More options are available to modify only if certain other options or variables satisfy certain conditions (_dependent options_).
_Dependent options_ can only be accessed if their precondition is satisfied. 

If the precondition is satisfied, they take on a default value, although they can be changed.
If the precondition is *not* satisfied, then their value is fixed and cannot be changed. For instance,

```bash
# in <top-level>/build
cmake .. -DSINGULARITY_USE_KOKKOS=OFF -DSINGULARITY_USE_CUDA=ON
```

Will have no effect (i.e. `SINGULARITY_USE_CUDA` will be set to `OFF`), because the precondition of `SINGULARITY_USE_CUDA` is for 
`SINGULARITY_USE_KOKKOS=ON`. 

Generally, _dependent options_ should only be used for specific use-cases where the defaults are not applicable.
For most scenarios, the preconditions and defaults are logically constructed and the most natural in practice 
(`SINGULARITY_TEST_*` are only available if `SINGLARITY_BUILD_TESTS` is enabled, for instance).

These options are listed in the following table, along with their preconditions:

|Option|Precondition|Default (condition true/false)| Comment |
|--|--|--|--|
| `SINGULARITY_USE_SPINER_WITH_HDF5` | `SINGULARITY_USE_SPINER=ON` | ON/OFF | Requests that `spiner` be configured for `HDF5` support.|
| `SINGULARITY_USE_CUDA`              | `SINGULARITY_USE_KOKKOS=ON`  | ON/OFF | Target nvidia GPUs for `Kokkos` offloading.|
| `SINGULARITY_USE_KOKKOSKERNELS`     | `SINGULARITY_USE_KOKKOS=ON`  | ON/OFF |Use Kokkos Kernels for linear algebra. Needed for mixed cell closure models on GPU.|
| `SINGULARITY_BUILD_CLOSURE` | `SINGULARITY_USE_KOKKOS=ON` `SINGULARITY_USE_KOKKOSKERNELS=ON` | ON/OFF | Mixed cell closure.|
| `SINGULARITY_BUILD_SESAME2SPINER` | `SINGULARITY_USE_SPINER=ON` `SINGULARITY_USE_SPINER_WITH_HDF5=ON` | ON/OFF | Builds the conversion tool sesame2spiner which makes files readable by SpinerEOS.|
| `SINGULARITY_BUILD_STELLARCOLLAPSE2SPINER` | `SINGULARITY_USE_SPINER=ON` `SINGULARITY_USE_SPINER_WITH_HDF5=ON` | ON/OFF | Builds the conversion tool stellarcollapse2spiner which optionally makes stellar collapse files faster to read.|
| `SINGULARITY_TEST_SESAME` | `SINGULARITY_BUILD_TESTS=ON` `SINGULARITY_BUILD_SESAME2SPINER=ON` | ON/OFF | Test the Sesame table readers.|
| `SINGULARITY_TEST_STELLAR_COLLAPSE` | `SINGULARITY_BUILD_TESTS=ON` `SINGULARITY_BUILD_STELLARCOLLAPSE2SPINER=ON`    | ON/OFF | Test the Stellar Collapse table readers.|
| `SINGULARITY_TEST_PYTHON`           | `SINGULARITY_BUILD_TESTS=ON` `SINGULARITY_BUILD_PYTHON=ON`    | ON/OFF    | Test the Python bindings.|

## CMake presets

To further aid the developer, `singularity-eos` is distributed with *Presets*, a list of common build options with naturally named 
labels that when used can reduce the need to input and remember the many options `singularity-eos` uses.
For a general overview of CMake presets, see the [cmake documentation on presets](https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html)

### Predefined presets 

Predefined presets are described with a `json` schema in the file `CMakePresets.json`. As an example:

```bash
# in <top-level>/build
$> cmake .. --preset="basic_with_testing"
Preset CMake variables:

  CMAKE_EXPORT_COMPILE_COMMANDS="ON"
  SINGULARITY_BUILD_TESTS="ON"
  SINGULARITY_USE_EOSPAC="ON"
  SINGULARITY_USE_SPINER="ON"

# ...
```

As you can see, CMake reports the configuration variables that the preset has used, and their values. A list of presets can be easily examined with:

```bash
# in <top-level>/build
$> cmake .. --list-presets
Available configure presets:

  "basic"
  "basic_with_testing"
  "kokkos_nogpu"
  "kokkos_nogpu_with_testing"
  "kokkos_gpu"
  "kokkos_gpu_with_testing"
```

When using presets, additional options may be readily appended to augment the required build.
For example, suppose that the `basic` preset is mostly sufficient, but you would like to enable building the closure models: 

```bash
# in <top-level>/build
$> cmake .. --preset="basic_with_testing" -DSINGULARITY_BUILD_CLOSURE=ON
# ...
```

### User defined presets 

The CMake preset functionality includes the ability of developers to define local presets in `CMakeUserPresets.json`.
`singularity-eos` explicitly does not track this file in Git, so developers can construct their own presets. 
All presets in the predefined `CMakePresets.json` are automatically included by CMake, so developers can 
build off of those if needed.

For instance, suppose you have a local checkout of the `kokkos` and `kokkos-kernels` codes that you're 
using to debug a GPU build, and you have these installed in `~/scratch/`.
Your `CMakeUserPresets.json` could look like:

```json
{
  "version": 1,
  "cmakeMinimumRequired": {
    "major": 3,
    "minor": 19
  },
  "configurePresets": [
    {
      "name": "my_local_build",
      "description": "submodule build using a local scratch install of kokkos",
      "inherits": [
        "kokkos_gpu_with_testing"
      ],
      "cacheVariables": {
        "Kokkos_DIR": "$env{HOME}/scratch/kokkos/lib/cmake/Kokkos",
        "KokkosKernels_DIR": "$env{HOME}/scratch/kokkoskernels/lib/cmake/KokkosKernels",
        "SINGULARITY_BUILD_PYTHON": "ON",
        "SINGULARITY_TEST_PYTHON": "OFF"
      }
    }
  ]
}
```

This inherits the predefined `kokkos_gpu_with_testing` preset, sets the `Kokkos*_DIR` cache variables to point `find_package()` 
to use these directories, and finally enables building the python bindings without including the python tests.

## Building in _submodule mode_

For _submodule mode_ to activate, a clone of the `singularity-eos` source should be placed
below the top-level of a host project

```bash
# An example directory layout when using singularity-eos in submodule mode
my_project
|_CMakeLists.txt
|_README.md 
|_src 
|_include
|_tpl/singularity-eos
```
`singularity-eos` is then imported using the `add_subdirectory()` command in CMake 

```cmake
# In your CMakeLists.txt
cmake_minimum_required(VERSION 3.19)
project(my_project)

add_executable(my_exec src/main.cc)
target_include_directories(my_exec include)

add_subdirectory(tpl/singularity-eos)

target_link_libraries(my_exec singularity-eos::singularity-eos)
```

This will expose the `singularity-eos` interface and library to your code, along with 
the interfaces of the internal dependencies

```c++
// in source of my_project 

#include<singularity-eos/eos/eos.hpp>
// from the internal ports-of-call submodule
#include<ports-of-call/portability>

// ...

using namespace singularity;
```

`singularity-eos` will build (along with internal dependencies) and be linked directly to your executable.

## Building in _standalone mode_ 

For _standalone_ mode, all required and optional dependencies are expected to be discoverable by CMake. This can be done several ways

1. (*preferred*) Use Spack to configure and install all the dependencies needed to build.
2. Use a system package manager (`apt-get`, `yum`, &t) to install dependencies.
3. Hand-build to a local filesystem, and configure your shell or CMake invocation to be aware of these installs

_standalone_ mode is the mode used to install `singularity-eos` to a system as a common library. If, for example, you use Spack to to install packages, `singularity-eos` will be built and installed in _standalone_ mode.



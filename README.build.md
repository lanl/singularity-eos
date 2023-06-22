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


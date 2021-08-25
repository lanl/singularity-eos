
#------------------------------------------------------------------------------#
# © 2021. Triad National Security, LLC. All rights reserved.  This
# program was produced under U.S. Government contract 89233218CNA000001
# for Los Alamos National Laboratory (LANL), which is operated by Triad
# National Security, LLC for the U.S.  Department of Energy/National
# Nuclear Security Administration. All rights in the program are
# reserved by Triad National Security, LLC, and the U.S. Department of
# Energy/National Nuclear Security Administration. The Government is
# granted for itself and others acting on its behalf a nonexclusive,
# paid-up, irrevocable worldwide license in this material to reproduce,
# prepare derivative works, distribute copies to the public, perform
# publicly and display publicly, and to permit others to do so.
#------------------------------------------------------------------------------#

#######################################
# ProcessConfiguration.cmake
#
# use provided configuration options
# to setup the build, including
# processing how dependencies are
# brought into the build.
#######################################

#######################################
## includes ###########################
#######################################

# Feature Summary used to present the user
# with an easy-to-read summary of options
# and features of the build configuration
#include(FeatureSummary)

# only present user options if preconditions are met
include(CMakeDependentOption)

# ${PROJECT_SOURCE_DIR}/cmake/ExternalConfig.cmake
# macros used to define the build configurations of
# submodules when selected
include(ExternalConfig)

#######################################
## options ############################
#######################################

# where available, use the internal source packages
# in ${PROJECT_SOURCE_DIR}/utils
option (SINGULARITY_USE_INTERNAL_DEPS "Use submodules when available" OFF)

# Options to build using Kokkos components
option (SINGULARITY_USE_KOKKOS "Use Kokkos for portability" OFF)
cmake_dependent_option (SINGULARITY_USE_CUDA "Enable cuda support" OFF "SINGULARITY_USE_KOKKOS" OFF)
cmake_dependent_option (SINGULARITY_USE_KOKKOSKERNELS "Use kokkos-kernels for linear algebra" ON "SINGULARITY_USE_KOKKOS" OFF)

# Options to build with EOSPAC
option (SINGULARITY_USE_EOSPAC "Pull in eospac" OFF)
cmake_dependent_option (SINGULARITY_INVERT_AT_SETUP "Use eospacs preinverted tables" OFF "SINGULARITY_USE_EOSPAC" OFF)

# Use HDF5 if available
option (SINGULARITY_USE_HDF5 "Pull in hdf5" ON)

# Build Fortran interface
option (SINGULARITY_USE_FORTRAN "Enable fortran bindings" ON)

# Build spiner components
option (SINGULARITY_BUILD_SESAME2SPINER "Compile sesame2spiner" OFF)
option (SINGULARITY_BUILD_STELLARCOLLAPSE2SPINER "Compile stellarcollapse2spiner" OFF)

# Options to control compiler output
option (SINGULARITY_BETTER_DEBUG_FLAGS
  "Better debug flags for singularity" ON)
option (SINGULARITY_HIDE_MORE_WARNINGS "hide more warnings" OFF)
option (SINGULARITY_BUILD_CLOSURE "Mixed cell closure" ON)

# Options for testing
option (SINGULARITY_BUILD_TESTS "Compile tests" OFF)
# TODO [mauneyc] Does sesame,stellar collapse require EOSPAC/spiner?
option (SINGULARITY_TEST_SESAME "Test the Sesame table readers" OFF)
option (SINGULARITY_TEST_STELLAR_COLLAPSE "Test the stellar collapse table readers" OFF)

# Misc and/or options I can't categorize
option (SINGULARITY_EOS_SKIP_EXTRAP "Turn off eospac extrap checks" OFF)
option (SINGULARITY_USE_SINGLE_LOGS "Use single precision logs. Can harm accuracy." OFF)
option (SINGULARITY_FMATH_USE_ORDER_4 "4th order interpolant for fast logs. Default is 7th order." OFF)
option (SINGULARITY_FMATH_USE_ORDER_5 "5th order interpolant for fast logs. Default is 7th order." OFF)

#######################################
## configuration logic ################
#######################################

singularity_message(STATUS "${ColorBold}Dependency setup${ColorReset}")
# Locate CUDA if requested
if(SINGULARITY_USE_CUDA)
    if(NOT TARGET CUDA::toolkit)
        find_package(CUDAToolkit REQUIRED)
    endif()
endif()
if(SINGULARITY_USE_EOSPAC)
    find_package(EOSPAC REQUIRED)
endif()
find_package(PortsofCall REQUIRED)

# The following dependencies may be taken from
# submodules (i.e. {top_dir}/utils/{dep}).
# see cmake/ExternalConfig.cmake for
# macro definition/use

#######################################
## Kokkos #############################
## and, if specified, KokkosKernels ###
#######################################
if(SINGULARITY_USE_KOKKOS)
    singularity_config_internal_kokkos(SINGULARITY_USE_CUDA)
    singularity_select_dep(
        DEP Kokkos
        GITURL https://github.com/kokkos/kokkos.git
        GITTAG origin/master
        SUBDIR kokkos
    )

    if(SINGULARITY_USE_KOKKOSKERNELS)
        singularity_config_internal_kokkoskernels()
        singularity_select_dep(
            DEP KokkosKernels
            GITURL https://github.com/kokkos/kokkos-kernels.git
            GITTAG origin/master
            SUBDIR ""
        )

    endif()
else()
    if(SINGULARITY_USE_CUDA)
        singularity_message(FATAL_ERROR "CUDA must be used with Kokkos")
    endif()
endif()

#######################################
## Eigen3 #############################
#######################################
if(NOT SINGULARITY_USE_KOKKOSKERNELS)
    singularity_select_dep(
        DEP Eigen3
        GITURL https://gitlab.com/libeigen/eigen.git
        GITTAG  origin/master
        SUBDIR "eigen"
    )
endif()

#######################################
## Catch2 #############################
#######################################
if(SINGULARITY_BUILD_TESTS)
    singularity_select_dep(
        DEP Catch2
        GITURL https://github.com/catchorg/Catch2.git
        GITTAG v2.13.4 # devel was broken
        SUBDIR "catch2"
    )
endif()

#######################################
## HDF5 ###############################
## NB:  CMake v3.21+ can correctly ####
##      handle `find_package(HDF5) ####
#######################################
if(SINGULARITY_USE_HDF5)
    find_package(HDF5 COMPONENTS C HL QUIET)

    if (HDF5_FOUND)
        add_library(${PROJECT_NAME}::hdf5 INTERFACE IMPORTED)
        set_target_properties(${PROJECT_NAME}::hdf5
        PROPERTIES
            INTERFACE_LINK_LIBRARIES "${HDF5_LIBRARIES};${HDF5_HL_LIBRARIES}"
            INTERFACE_COMPILE_DEFINITIONS "SINGULARITY_USE_HDF5"
            INTERFACE_INCLUDE_DIRECTORIES "${HDF5_INCLUDE_DIRS}"
        )
        if(HDF5_IS_PARALLEL)
            # TODO: this shouldn't be necessary, tho
            # i'm leaving it in for now
            if(NOT TARGET MPI::MPI_CXX)
                find_package(MPI COMPONENTS CXX QUIET)
                if(NOT MPI_FOUND)
                    message(FATAL_ERROR "The HDF5 package requires MPI")
                endif()
            endif()
        endif()
    endif()
endif()

# mauneyc: is this only for tests?
if (SINGULARITY_USE_HDF5 AND SINGULARITY_BUILD_STELLARCOLLAPSE2SPINER)
   add_subdirectory(${PROJECT_SOURCE_DIR}/stellarcollapse2spiner)
   install(TARGETS stellarcollapse2spiner DESTINATION bin)
endif()

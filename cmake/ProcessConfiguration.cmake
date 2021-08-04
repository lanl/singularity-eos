# Feature Summary used to present the user
# with an easy-to-read summary of options
# and features of the build configuration
include(FeatureSummary)

# options
option (SINGULARITY_USE_SUBMODULES "Default to submodule dependencies" OFF)
option (SINGULARITY_USE_HDF5 "Pull in hdf5" ON)
option (SINGULARITY_USE_FORTRAN "Enable fortran bindings" ON)
option (SINGULARITY_USE_KOKKOS "Use Kokkos for portability" OFF)
option (SINGULARITY_USE_EOSPAC "Pull in eospac" OFF)
option (SINGULARITY_USE_CUDA "Enable cuda support" OFF)
option (SINGULARITY_USE_KOKKOSKERNELS
  "Use kokkos-kernels for linear algebra" OFF)
option (SINGULARITY_INVERT_AT_SETUP "Use eospacs preinverted tables" OFF)
option (SINGULARITY_BUILD_TESTS "Compile tests" OFF)
option (SINGULARITY_BUILD_SESAME2SPINER "Compile sesame2spiner" OFF)
option (SINGULARITY_BUILD_STELLARCOLLAPSE2SPINER "Compile stellarcollapse2spiner" OFF)
option (SINGULARITY_BETTER_DEBUG_FLAGS
  "Better debug flags for singularity" ON)
option (SINGULARITY_HIDE_MORE_WARNINGS "hide more warnings" OFF)
option (SINGULARITY_SUBMODULE_MODE "Submodule mode" OFF)
option (SINGULARITY_BUILD_CLOSURE "Mixed cell closure" ON)
option (SINGULARITY_TEST_SESAME "Test the Sesame table readers" OFF)
option (SINGULARITY_TEST_STELLAR_COLLAPSE "Test the stellar collapse table readers" OFF)

if(SINGULARITY_USE_CUDA)
    if(NOT TARGET CUDA::toolkit)
        find_package(CUDAToolkit REQUIRED)
    endif()
endif()

if(SINGULARITY_USE_SUBMODULES)
    include(SubmoduleConfig)
endif()

if(SINGULARITY_USE_KOKKOS)
    if(NOT TARGET Kokkos::kokkos)
        find_package(Kokkos REQUIRED)
    endif()
else()
    if(SINGULARITY_USE_CUDA)
        message(FATAL_ERROR "CUDA must be used with Kokkos")
    endif()
    if(SINGULARITY_USE_KOKKOSKERNELS)
        message(FATAL_ERROR "KokkosKernels must be used with Kokkos")
    endif()
endif()

if(SINGULARITY_USE_KOKKOSKERNELS)
    if(NOT Kokkos::kokkoskernels)
        find_package(KokkosKernels REQUIRED)
    endif()
else()
    if(NOT Eigen3::Eigen)
        find_package(Eigen3 REQUIRED)
    endif()
endif()

if(SINGULARITY_USE_EOSPAC)
    find_package(EOSPAC)
endif()

find_package(PortsofCall REQUIRED)
find_package(Catch2 REQUIRED)

# Setup HDF5
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
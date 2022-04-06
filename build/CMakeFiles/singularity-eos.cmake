# config options


set(SINGULARITY_BUILD_SHARED_LIBS "ON")
set(SINGULARITY_EOS_SKIP_EXTRAP "OFF")
set(SINGULARITY_INVERT_AT_SETUP "OFF")
set(SINGULARITY_USE_CUDA "OFF")
set(SINGULARITY_USE_EIGEN "OFF")
set(SINGULARITY_USE_EOSPAC "OFF")
set(SINGULARITY_USE_FMATH_ORDER "7")
set(SINGULARITY_USE_FMATH_ORDERS "7;5;0")
set(SINGULARITY_USE_HDF5 "OFF")
set(SINGULARITY_USE_KOKKOS "OFF")
set(SINGULARITY_USE_KOKKOSKERNELS "OFF")
set(SINGULARITY_USE_SINGLE_LOGS "OFF")

if(NOT TARGET singularity-eos AND NOT singularity-eos_BINARY_DIR)
  include(/usr/local/lib/cmake/singularity-eos/singularity-eosTargets.cmake)
endif()

# handle TPLs

find_package(singularity-eosCMake REQUIRED)
list(APPEND CMAKE_MODULE_PATH ${singularity-eos_PREFIX})

set(SINGULARITY_USE_KOKKOS OFF)
include(kokkos)

set(SINGULARITY_USE_KOKKOSKERNELS OFF)
include(kokkoskernels)

set(SINGULARITY_USE_CUDA OFF)
include(cuda)

set(SINGULARITY_USE_EOSPAC) OFF)
include(eospac)

set(SINGULARITY_USE_HDF5) OFF)
include(hdf5)

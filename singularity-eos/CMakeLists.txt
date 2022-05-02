#------------------------------------------------------------------------------#
# © 2021-2022. Triad National Security, LLC. All rights reserved.  This
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

set(EOS_HEADERS
    base/constants.hpp
    base/eos_error.hpp
    base/fast-math/logs.hpp
    base/root-finding-1d/root_finding.hpp
    eos/eos.hpp
    eos/eos_builder.hpp
    eos/eos_variant.hpp
    eos/modifiers/eos_unitsystem.hpp
    eos/modifiers/scaled_eos.hpp
    eos/modifiers/shifted_eos.hpp
    eos/modifiers/relativistic_eos.hpp
)

set(EOS_SRCS
<<<<<<< HEAD:singularity-eos/CMakeLists.txt
    eos/eos_spiner.cpp
    eos/eos_jwl.cpp
    eos/eos_builder.cpp
    eos/eos_ideal.cpp
    eos/eos_gruneisen.cpp
    eos/eos_davis.cpp
    eos/eos_stellar_collapse.cpp
    eos/eos_eospac.cpp
)

if (SINGULARITY_BUILD_CLOSURE)
  list(APPEND EOS_SRCS eos/singularity_eos.cpp)
  list(APPEND EOS_HEADERS eos/singularity_eos.hpp closure/mixed_cell_models.hpp)
endif()

if (SINGULARITY_USE_FORTRAN)
  list(APPEND EOS_SRCS eos/singularity_eos.f90)
endif()

# export to parent scope
set(EOS_HEADERS ${EOS_HEADERS} PARENT_SCOPE)
set(EOS_SRCS ${EOS_SRCS} PARENT_SCOPE)

#------------------------------------------------------------------------------#
# © 2021 - 2022. Triad National Security, LLC. All rights reserved.  This
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

if(SINGULARITY_USE_FORTRAN AND SINGULARITY_BUILD_CLOSURE)
  add_library(eos_f OBJECT singularity_eos.f90)
endif()

set(EOS_SRCS
    eos.hpp
    modifiers/eos_unitsystem.hpp
    modifiers/scaled_eos.hpp
    modifiers/shifted_eos.hpp
    modifiers/relativistic_eos.hpp
    eos_spiner.cpp
    eos_jwl.cpp
    eos_builder.cpp
    eos_ideal.cpp
    eos_gruneisen.cpp
    eos_davis.cpp
    eos_stellar_collapse.cpp
    eos_eospac.cpp
    )
if (SINGULARITY_BUILD_CLOSURE)
  list(APPEND EOS_SRCS singularity_eos.cpp)
endif()

add_library(eos_objs OBJECT ${EOS_SRCS})

target_include_directories(eos_objs PUBLIC
            ${CMAKE_CURRENT_BINARY_DIR}
            ${CMAKE_CURRENT_SOURCE_DIR}
            )

target_link_libraries(eos_objs singularity-eos::libs singularity-eos::flags)

if(SINGULARITY_USE_FORTRAN AND SINGULARITY_BUILD_CLOSURE)
  add_library(eos $<TARGET_OBJECTS:eos_objs> $<TARGET_OBJECTS:eos_f>)
else()
  add_library(eos $<TARGET_OBJECTS:eos_objs>)
endif()

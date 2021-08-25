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

project(eos)

add_library(${PROJECT_NAME})
add_library(singularity-eos::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

target_sources(${PROJECT_NAME}
PRIVATE
    eos_spiner.cpp
    eos_jwl.cpp
    eos_builder.cpp
    eos_ideal.cpp
    eos_gruneisen.cpp
    eos_davis.cpp
    eos_stellar_collapse.cpp
    eos_eospac.cpp
    $<$<BOOL:${SINGULARITY_BUILD_CLOSURE}>:singularity_eos.cpp>
    $<$<BOOL:${SINGULARITY_USE_FORTRAN}>:singularity_eos.f90>
)

# TODO: double-check INSTALL_INTERFACE
target_include_directories(${PROJECT_NAME}
PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)

target_link_libraries(${PROJECT_NAME} PUBLIC singularity-eos::flags)

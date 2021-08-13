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



add_library(eos)
add_library(singularity-eos::eos ALIAS eos)

target_sources(eos
PRIVATE
    src/eos_spiner.cpp
    src/eos_jwl.cpp
    src/eos_builder.cpp
    src/eos_ideal.cpp
    src/eos_gruneisen.cpp
    src/eos_davis.cpp
    src/eos_stellar_collapse.cpp
    src/eos_eospac.cpp
    $<$<BOOL:${SINGULARITY_BUILD_CLOSURE}>:src/singularity_eos.cpp>
    $<$<BOOL:${SINGULARITY_USE_FORTRAN}>:src/singularity_eos.f90>
#PUBLIC
#    include/eos/eos.hpp
#    include/eos/eos_builder.hpp
#    include/eos/eos_variant.hpp
#    include/eos/singularity_eos.hpp
#    include/eos/base/constants.hpp
#    include/eos/base/eos_error.hpp
#    include/eos/closure/mixed_cell_models.hpp
)

target_include_directories(eos
PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

target_link_libraries(eos ${PROJECT_NAME}::flags)

set_target_properties(eos PROPERTIES EXPORT_NAME singularity-eos) 

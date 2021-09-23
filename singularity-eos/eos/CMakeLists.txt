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

set(eos_HEADERS
  eos_builder.hpp
  eos_variant.hpp
  eos.hpp
  singularity_eos.hpp
  base/constants.hpp
  base/eos_error.hpp
  modifiers/eos_unitsystem.hpp
  modifiers/scaled_eos.hpp
  modifiers/shifted_eos.hpp
  modifiers/relativistic_eos.hpp
)

set(eos_SOURCES
  eos_spiner.cpp
  eos_jwl.cpp
  eos_builder.cpp
  eos_ideal.cpp
  eos_gruneisen.cpp
  eos_davis.cpp
  eos_stellar_collapse.cpp
  eos_eospac.cpp
)
if(SINGULARITY_BUILD_CLOSURE)
  set(eos_SOURCES
    ${eos_SOURCES}
    singularity_eos.cpp
  )
endif()
if(SINGULARITY_USE_FORTRAN)
  set(eos_SOURCES
    ${eos_SOURCES}
    singularity_eos.f90
  )

endif()

#------------------------------------------------------------------------------#
# Export header list to parent scope.
#------------------------------------------------------------------------------#

set(eos_HEADERS
  ${eos_HEADERS}
  PARENT_SCOPE
)

#------------------------------------------------------------------------------#
# Export source list to parent scope.
#------------------------------------------------------------------------------#

set(eos_SOURCES
  ${eos_SOURCES}
  PARENT_SCOPE
)
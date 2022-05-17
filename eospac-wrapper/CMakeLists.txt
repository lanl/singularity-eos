#------------------------------------------------------------------------------
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
#------------------------------------------------------------------------------

option(SINGULARITY_EOSPAC_SKIP_EXTRAP "Skip extrapolation checks" OFF)

add_library(eospac-wrapper
  eospac_wrapper.hpp
  eospac_wrapper.cpp
)

if (SINGULARITY_EOSPAC_SKIP_EXTRAP)
  target_compile_definitions(eospac-wrapper
    PRIVATE
      SINGULARITY_EOSPAC_SKIP_EXTRAP
  )
endif() #SINGULARITY_EOSPAC_SKIP_EXTRAP

target_link_libraries(eospac-wrapper 
  PUBLIC
    EOSPAC::eospac
)

#TODO: should this be here? not sure what's appropriate 
# (could put line in top-level CMakeLists.txt)
#install(TARGETS eospac-wrapper EXPORT singularity-eosTargets DESTINATION lib)
list(APPEND SINGULARITY_EXPORT_TARGETS eospac-wrapper)

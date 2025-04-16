#------------------------------------------------------------------------------
# © 2021-2023. Triad National Security, LLC. All rights reserved.  This
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

add_executable(get_sound_speed_press
  get_sound_speed_press.cpp)
target_link_libraries(get_sound_speed_press PRIVATE
  singularity-eos::singularity-eos)

if (SINGULARITY_USE_SPINER_WITH_HDF5)
  add_executable(eos_grid eos_grid.cpp)
  target_link_libraries(eos_grid PRIVATE
    singularity-eos::singularity-eos)
endif()

if(SINGULARITY_USE_EOSPAC AND SINGULARITY_USE_SPINER_WITH_HDF5)
  add_executable(get_sesame_state
    get_sesame_state.cpp)
  target_link_libraries(get_sesame_state PRIVATE
    singularity-eos::singularity-eos)
  add_executable(get_spiner_bounds
    get_spiner_bounds.cpp)
  target_link_libraries(get_spiner_bounds PRIVATE
    singularity-eos::singularity-eos)
endif()

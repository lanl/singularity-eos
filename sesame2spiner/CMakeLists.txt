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

add_library(sesame2spiner-lib
  io_eospac.cpp
  io_eospac.hpp
  
  generate_files.cpp
  generate_files.hpp
  
  parse_cli.cpp
  parse_cli.hpp
  
  parser.cpp
  parser.hpp
  )

target_include_directories(sesame2spiner-lib PUBLIC  
                           ${CMAKE_CURRENT_BINARY_DIR}
                           ${CMAKE_CURRENT_SOURCE_DIR}
                          )

target_link_libraries(sesame2spiner-lib
  singularity-eos::libs singularity-eos::flags)

add_executable(sesame2spiner
  main.cpp
  )

target_link_libraries(sesame2spiner
  sesame2spiner-lib
  singularity-eos::libs singularity-eos::flags)

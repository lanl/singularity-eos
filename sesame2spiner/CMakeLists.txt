#------------------------------------------------------------------------------
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
#------------------------------------------------------------------------------

add_executable(sesame2spiner
               io_eospac.cpp
               #test.cpp
               generate_files.cpp
               parser.cpp
               main.cpp
              )

target_include_directories(sesame2spiner PUBLIC  
                           ${CMAKE_CURRENT_BINARY_DIR}
                           ${CMAKE_CURRENT_SOURCE_DIR}
                          )

target_link_libraries(sesame2spiner nlohmann_json::nlohmann_json singularity-eos::libs singularity-eos::flags)

# TODO: Add tests for sesame2spiner here.

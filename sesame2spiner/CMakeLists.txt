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
project(sesame2spiner)
find_package(singularity-eos REQUIRED)

add_executable(${PROJECT_NAME})
target_sources(${PROJECT_NAME}
PRIVATE
    io_eospac.cpp
	io_eospac.hpp
    generate_files.cpp
	generate_files.hpp
    parse_cli.cpp
    parse_cli.hpp
    parser.cpp
    parser.hpp
PUBLIC
    main.cpp
)

target_include_directories(${PROJECT_NAME} PUBLIC
                            "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
                            "$<INSTALL_INTERFACE:${CMAKE_INSTALL_BINDIR}>"
                          )

target_link_libraries(${PROJECT_NAME} PUBLIC singularity-eos::singularity-eos)
install(TARGETS ${PROJECT_NAME} DESTINATION bin)
# TODO: Add tests for sesame2spiner here.

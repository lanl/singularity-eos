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

target_include_directories(${PROJECT_NAME}
PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)

target_link_libraries(${PROJECT_NAME} PUBLIC singularity-eos::flags)

#find_package(Python COMPONENTS Interpreter Development)
#find_package(pybind11 CONFIG)

add_subdirectory(${PROJECT_SOURCE_DIR}/utils/pybind11 pybind11)

find_package(Python REQUIRED COMPONENTS Development)
pybind11_add_module(singularity_eos module.cpp)
target_link_libraries(singularity_eos PRIVATE singularity-eos singularity-eos::flags singularity-eos::libs)

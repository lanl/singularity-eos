project(closure)

add_library(${PROJECT_NAME} INTERFACE)
add_library(singularity-eos::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

target_include_directories(
  ${PROJECT_NAME}
  INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)


#target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}::flags)
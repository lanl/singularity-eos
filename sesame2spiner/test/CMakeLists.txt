# Native-table tests exercise the conversion library independently of EOS models.
add_executable(sesame2spiner_native_tests
  ${PROJECT_SOURCE_DIR}/test/catch2_define.cpp
  test_native_tables.cpp)
target_link_libraries(sesame2spiner_native_tests PRIVATE
  Catch2::Catch2 sesame2spiner::sesame2spiner)
include(Catch)
catch_discover_tests(sesame2spiner_native_tests PROPERTIES TIMEOUT 120)

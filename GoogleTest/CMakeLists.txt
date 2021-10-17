add_executable(unit_test test.cpp)
find_package(GTest REQUIRED)
target_compile_features(unit_test PUBLIC cxx_std_17)
target_link_libraries(unit_test GTest::GTest GTest::Main)

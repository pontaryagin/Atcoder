add_executable(unit_test test.cpp)
find_package(GTest REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_DEBUG")
target_compile_features(unit_test PUBLIC cxx_std_17)
target_link_libraries(unit_test GTest::GTest GTest::Main)
target_include_directories(unit_test
  PUBLIC ${PROJECT_SOURCE_DIR}/ac-library ${PROJECT_SOURCE_DIR}/src
  )

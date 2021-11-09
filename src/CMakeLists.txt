add_executable(main.out main.cpp)
target_compile_features(main.out PUBLIC cxx_std_17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_DEBUG")
target_include_directories(main.out
  PUBLIC ${PROJECT_SOURCE_DIR}/ac-library
  )

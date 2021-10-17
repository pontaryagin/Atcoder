add_executable(a.out main.cpp)
target_compile_features(a.out PUBLIC cxx_std_17)
target_include_directories(a.out
  PUBLIC ${PROJECT_SOURCE_DIR}/atcoder
  )

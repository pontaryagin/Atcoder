add_executable(main.out main.cpp)
target_compile_features(main.out PUBLIC cxx_std_17)
target_include_directories(main.out
  PUBLIC ${PROJECT_SOURCE_DIR}/ac-library
  )

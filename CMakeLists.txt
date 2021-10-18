cmake_minimum_required(VERSION 3.16)
project(AtCoder CXX)
# サブディレクトリを登録
add_subdirectory(Project1)
add_subdirectory(GoogleTest)

# set(CMAKE_SOURCE_DIR ${CMAKE_HOME_DIRECTORY})
# set(CMAKE_BINARY_DIR ${CMAKE_HOME_DIRECTORY}/build)
set(DEFAULT_OUT_OF_SOURCE_FOLDER "build")

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/debug)

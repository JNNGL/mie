cmake_minimum_required(VERSION 3.10)

project(libmie)

set(CMAKE_CXX_STANDARD 17)

file(GLOB_RECURSE LIBMIE_SOURCES *.cpp)
add_library(libmie STATIC ${LIBMIE_SOURCES})
target_include_directories(libmie PUBLIC ${PROJECT_SOURCE_DIR})
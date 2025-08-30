cmake_minimum_required(VERSION 3.10)

project(bakemie)

set(CMAKE_CXX_STANDARD 17)

file(GLOB_RECURSE BAKEMIE_SOURCES *.cpp)
add_executable(bakemie ${BAKEMIE_SOURCES})
target_include_directories(bakemie PUBLIC ${PROJECT_SOURCE_DIR})
target_link_libraries(bakemie PUBLIC libmie)
cmake_minimum_required(VERSION 3.30)
project(lorenz_mie)

set(CMAKE_CXX_STANDARD 17)

include_directories(${CMAKE_SOURCE_DIR})

add_executable(lorenz_mie main.cpp lib/mie.cpp)
target_link_libraries(lorenz_mie matplot)

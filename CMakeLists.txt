cmake_minimum_required(VERSION 3.30)
project(lorenz_mie LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)

include_directories(${CMAKE_SOURCE_DIR})

set(MATPLOTPP_BUILD_EXPERIMENTAL_OPENGL_BACKEND ON)
add_subdirectory(thirdparty/matplot)

option(ENABLE_GPU "Enable GPU acceleration" ON)

set(SOURCE_FILES main.cpp lib/mie.cpp)

add_executable(lorenz_mie ${SOURCE_FILES})
target_link_libraries(lorenz_mie matplot_opengl)

if (ENABLE_GPU)
    target_compile_definitions(lorenz_mie PUBLIC ENABLE_GPU)
    set(TARGET_LANGUAGE CUDA)
    message("-- Using GPU acceleration")
else ()
    set(TARGET_LANGUAGE CXX)
endif ()

set_source_files_properties(${SOURCE_FILES} PROPERTIES LANGUAGE ${TARGET_LANGUAGE})
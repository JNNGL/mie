cmake_minimum_required(VERSION 3.10)

project(libmie)

include(CheckLanguage)

option(ENABLE_CUDA_BACKEND "Enable CUDA backend" ON)

set(CMAKE_CXX_STANDARD 17)

file(GLOB_RECURSE LIBMIE_SOURCES *.cpp)
add_library(libmie STATIC ${LIBMIE_SOURCES})
target_include_directories(libmie PUBLIC ${PROJECT_SOURCE_DIR})

if (ENABLE_CUDA_BACKEND)
    check_language(CUDA)

    if (CMAKE_CUDA_COMPILER)
        enable_language(CUDA)

        file(GLOB_RECURSE LIBMIE_CUDA_SOURCES *.cu)
        target_sources(libmie PUBLIC ${LIBMIE_CUDA_SOURCES})
        target_compile_definitions(libmie PUBLIC ENABLE_CUDA_BACKEND)

        message(STATUS "Enabled CUDA backend")
    else()
        message(WARNING "Disabled CUDA backend: couldn't find CUDA compiler")
    endif()
endif()
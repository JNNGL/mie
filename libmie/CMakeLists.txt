cmake_minimum_required(VERSION 3.18)

project(libmie)

option(ENABLE_OPENCL_BACKEND "Enable OpenCL backend" ON)

set(CMAKE_CXX_STANDARD 17)

file(GLOB_RECURSE LIBMIE_SOURCES *.cpp)
add_library(mie STATIC ${LIBMIE_SOURCES})
target_include_directories(mie PUBLIC ${PROJECT_SOURCE_DIR})

if (ENABLE_OPENCL_BACKEND)
    find_package(OpenCL)

    if (OpenCL_FOUND)
        target_link_libraries(mie PRIVATE OpenCL::OpenCL)
        target_compile_definitions(mie PUBLIC LIBMIE_OPENCL_BACKEND)

        set(EMBED_HEADERS_LOCATION "${CMAKE_BINARY_DIR}/embed")

        file(READ "${PROJECT_SOURCE_DIR}/mie/backend/builtin/opencl_kernel.c" BYTES HEX)
        string(REGEX REPLACE "(..)" "0x\\1, " BYTES "${BYTES}")

        string(CONFIGURE "const std::uint8_t opencl_kernel_data[] = { ${BYTES} 0 };" DATA_DECLARATION)

        set(EMBED_HEADER "\
            #pragma once

            #include <cstdint>

            ${DATA_DECLARATION}
        ")

        string(CONFIGURE "${BASIC_HEADER}" BASIC_HEADER)

        file(WRITE "${EMBED_HEADERS_LOCATION}/opencl_kernel.h" "${EMBED_HEADER}")
        target_include_directories(mie PRIVATE "${EMBED_HEADERS_LOCATION}")
        target_compile_definitions(mie PRIVATE CL_TARGET_OPENCL_VERSION=200)

        message(STATUS "Enabled OpenCL backend")

        add_library(opencl_kernel_dummy STATIC "${PROJECT_SOURCE_DIR}/mie/backend/builtin/opencl_kernel.c")
        target_compile_definitions(opencl_kernel_dummy PRIVATE cl_khr_fp64 CMAKE_PROJECT)
    else()
        message(STATUS "Disabled OpenCL backend: couldn't find OpenCL")
    endif()
endif()

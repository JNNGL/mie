cmake_minimum_required(VERSION 3.30)
project(lorenz_mie)

option(ENABLE_GPU "Enable GPU acceleration" ON)
option(BUILD_PLOT "Build `plot_mie` executable" ON)

if (ENABLE_GPU)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
endif()

set(CMAKE_CXX_STANDARD 17)

include_directories(${CMAKE_SOURCE_DIR})

add_subdirectory(thirdparty/matplot)

set(LIB_SOURCE_FILES lib/mie.cpp)
set(PLOT_MIE_SOURCE_FILES plot.cpp)

add_library(mie STATIC ${LIB_SOURCE_FILES})

add_executable(plot_mie ${PLOT_MIE_SOURCE_FILES})
target_link_libraries(plot_mie PUBLIC matplot mie)

if (ENABLE_GPU)
    add_compile_definitions(ENABLE_GPU)
    set(TARGET_LANGUAGE CUDA)
    message("-- Using GPU acceleration")
else ()
    set(TARGET_LANGUAGE CXX)
endif ()

set_source_files_properties(${PLOT_MIE_SOURCE_FILES} ${LIB_SOURCE_FILES} PROPERTIES LANGUAGE ${TARGET_LANGUAGE})
cmake_minimum_required(VERSION 3.18)

project(bakemie)

set(CMAKE_CXX_STANDARD 17)

file(GLOB_RECURSE BAKEMIE_SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp)
add_executable(bakemie ${BAKEMIE_SOURCES})
target_include_directories(bakemie PUBLIC ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(bakemie PUBLIC mie)

configure_file(${PROJECT_SOURCE_DIR}/resources/zap-vga09.psf ${PROJECT_BINARY_DIR}/font.psf COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/resources/RobotoMono-Regular.ttf ${PROJECT_BINARY_DIR}/font.ttf COPYONLY)

if(EXISTS ${PROJECT_SOURCE_DIR}/third_party/stb AND IS_DIRECTORY ${PROJECT_SOURCE_DIR}/third_party/stb)
    target_include_directories(bakemie PUBLIC ${PROJECT_SOURCE_DIR}/third_party/stb)
    target_compile_definitions(bakemie PUBLIC HAS_STB)
endif()
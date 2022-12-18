cmake_minimum_required(VERSION 3.16)
project(CS171-assignment4 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)

if(CMAKE_BUILD_TYPE STREQUAL Debug)
    add_compile_definitions(MY_DEBUG)
endif()

find_package(OpenMP REQUIRED)

add_subdirectory(libs)
add_subdirectory(src)

add_executable(${PROJECT_NAME}-main main.cpp)

target_link_libraries(${PROJECT_NAME}-main
        PRIVATE
        renderer)
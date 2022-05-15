cmake_minimum_required(VERSION 3.18)

project(CTFD1)

# gnuplot is required for visualization
find_program(GNUPLOT NAMES gnuplot REQUIRED True)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)

add_subdirectory(chapter_03)

add_executable(gp_test gp_test.c)

target_link_libraries(gp_test m)

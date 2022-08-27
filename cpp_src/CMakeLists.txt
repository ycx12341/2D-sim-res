cmake_minimum_required(VERSION 3.23)
project(sim_2d_cpp VERSION 0.0.1)

set(CMAKE_CXX_STANDARD 17)

include_directories(eigen-3.4.0)

add_executable(cpp_src main.cpp algo/calculate_sse.h algo/calculate_sse.cpp algo/generate_pattern.cpp seed/seed.h seed/seed.cpp seed/sample.cpp config.h collection/collection.h collection/matrix_ext.h algo/solve_pde.cpp)
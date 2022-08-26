cmake_minimum_required(VERSION 3.23)
project(c_src C)

set(CMAKE_C_STANDARD 11)

include_directories(util/math)

add_executable(c_src
        config.h
        main.c
        seed/seed.h
        config.h seed/seed.c seed/sample.c collection/collection.h collection/arraylist.c collection/collection.c algo/algo.h algo/calculate_scc.c algo/generate_pattern.c algo/pars.c matrix/matrix.h matrix/matrix.c algo/proliferation.c algo/pde.c)

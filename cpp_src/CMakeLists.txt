cmake_minimum_required(VERSION 3.23)
project(sim_2d_cpp VERSION 0.0.1)

set(CMAKE_CXX_STANDARD 17)

add_executable(
        cpp_src
        main.cpp
        algo/calculate.h
        algo/pattern.h
        seed/seed.h
        seed/seed.cpp
        seed/sample.cpp
        config.h
        collection/collection.h
        algo/pde.h
        collection/matrix.h
        algo/scc.h
        algo/sim2d.h
        seed/nrom.cpp
        algo/pars.h
        ref_den/t3_ref_den.h
        sample/sample.h
        sample/scc.cpp
)
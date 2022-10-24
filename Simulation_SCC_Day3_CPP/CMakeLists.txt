cmake_minimum_required(VERSION 3.23)
project(sim_2d_cpp VERSION 0.0.1)

set(CMAKE_CXX_STANDARD 17)

add_executable(
        Simulation_SCC_Day3_CPP
        main.cpp
        config.h

        algo/calculate.h
        algo/pars.h
        algo/pattern.h
        algo/pde.h
        algo/scc.h
        algo/sim2d.h

        collection/collection.h
        collection/matrix.h

        ref_den/ref_den.h

        sample/sample.h
        sample/scc.cpp

        seed/nrom.cpp
        seed/sample.cpp
        seed/seed.cpp
        seed/seed.h
)
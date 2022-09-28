cmake_minimum_required(VERSION 3.23)
project(SciHPC)

set(CMAKE_CXX_STANDARD 23)
add_executable(SciHPC main.cpp scalar_data.cpp scalar_data.h global.h structured_grid.cpp structured_grid.h
        matrix_solver.cpp matrix_solver.h vector_data.cpp vector_data.h global.cpp global.h runge_kutta.cpp runge_kutta.h boundary_condition.cpp boundary_condition.h)

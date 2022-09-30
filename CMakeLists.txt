cmake_minimum_required(VERSION 3.23)
project(SciHPC)

set(CMAKE_CXX_STANDARD 23)
add_executable(SciHPC main.cpp scihpc/scalar_data.cpp scihpc/scalar_data.h scihpc/global.h scihpc/structured_grid.cpp scihpc/structured_grid.h
        scihpc/matrix_solver.cpp scihpc/matrix_solver.h scihpc/vector_data.cpp scihpc/vector_data.h scihpc/global.cpp scihpc/global.h scihpc/runge_kutta.cpp scihpc/runge_kutta.h scihpc/boundary_condition.cpp scihpc/boundary_condition.h scihpc/ccd.cpp scihpc/ccd.h)

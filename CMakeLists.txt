cmake_minimum_required(VERSION 3.23)
project(SciHPC)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
set(CMAKE_CXX_STANDARD 23)

include_directories(scihpc)

add_executable(SciHPC
#        main.cpp
#        benchmarks_convergence.cpp
        benchmarks_inviscid_Burger's_equation.cpp
#        benchmark_2d_vortex.cpp
        scihpc/boundary_condition.cpp
        scihpc/boundary_condition.h
        scihpc/flux.cpp
        scihpc/flux.h
        scihpc/global.cpp
        scihpc/global.h
        scihpc/lsm.cpp
        scihpc/lsm.h
        scihpc/matrix_solver.cpp
        scihpc/matrix_solver.h
        scihpc/runge_kutta.cpp
        scihpc/runge_kutta.h
        scihpc/scalar_data.cpp
        scihpc/scalar_data.h
        scihpc/source.cpp
        scihpc/source.h
        scihpc/structured_grid.cpp
        scihpc/structured_grid.h
        scihpc/vector_data.cpp
        scihpc/vector_data.h
        scihpc/vtkWriter.cpp
        scihpc/vtkWriter.h
        scihpc/ccd_base.cpp
        scihpc/ccd_base.h
        scihpc/ccd_solver.cpp
        scihpc/ccd_solver.h
        scihpc/uccd_base.cpp
        scihpc/uccd_base.h
        scihpc/uccd_solver.cpp
        scihpc/uccd_solver.h
        scihpc/derivatives_solver.h
        scihpc/wrapper.cpp
        scihpc/wrapper.h scihpc/dummy_data.h scihpc/dummy_data.cpp scihpc/weno_solver.cpp scihpc/weno_solver.h scihpc/derivatives_solver.cpp)

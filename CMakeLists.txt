cmake_minimum_required(VERSION 3.23)
project(SciHPC)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
set(CMAKE_CXX_STANDARD 23)

include_directories(scihpc)

add_executable(SciHPC
#        main.cpp
        benchmarks_1d_wave.cpp
#        benchmarks_1d_IBE.cpp
#        benchmarks_2d_vortex.cpp
#        benchmarks_2d_zalesak.cpp
#        benchmarks_lsf_initialization.cpp
#        benchmarks_lsf_curvature.cpp.cp
#        benchmarks_lid_driven_cavity.cpp
        scihpc/simple_bc.cpp
        scihpc/simple_bc.h
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
        scihpc/wrapper.h
        scihpc/dummy_data.h
        scihpc/dummy_data.cpp
        scihpc/weno_solver.cpp
        scihpc/weno_solver.h
        scihpc/derivatives_solver.cpp
        scihpc/godunov_gradient.cpp
        scihpc/godunov_gradient.h
        scihpc/wrapper_func.cpp
        scihpc/wrapper_func.h
        scihpc/projection_method.cpp
        scihpc/projection_method.h
        scihpc/no_slip.cpp
        scihpc/no_slip.h scihpc/velocities_bc.cpp scihpc/velocities_bc.h scihpc/slip.cpp scihpc/slip.h)

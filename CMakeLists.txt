cmake_minimum_required(VERSION 3.21.0)
project(SciHPC)

if(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp -O3")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O3")
endif()

set(CMAKE_CXX_STANDARD 17)

include_directories(scihpc)

add_executable(SciHPC
#        main.cpp
#        benchmarks_1d_wave.cpp
#        benchmarks_1d_IBE.cpp
#        benchmarks_2d_vortex.cpp
#        benchmarks_2d_zalesak.cpp
#        benchmarks_lsf_initialization.cpp
#        benchmarks_lsf_curvature.cpp
#        benchmarks_2d_cavity_flows.cpp
#        benchmarks_2d_dambreak.cpp
#        benchmarks_2d_RayleighTaylorInstability.cpp
#        benchmarks_2d_bubble_rising.cpp
#        benchmarks_2d_multigrid_Poisson.cpp
#        benchmarks_3d_dambreak.cpp
#        benchmarks_3d_multigrid_Poisson.cpp
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
        scihpc/SharedSolvers.h
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
        scihpc/no_slip.h
        scihpc/bc_factory.cpp
        scihpc/bc_factory.h
        scihpc/slip.cpp
        scihpc/slip.h
        scihpc/dirichlet.cpp
        scihpc/dirichlet.h
        scihpc/neumann.cpp
        scihpc/neumann.h
        scihpc/periodic.cpp
        scihpc/periodic.h
        scihpc/second_order_solver.cpp
        scihpc/second_order_solver.h
        scihpc/multigrid.cpp
        scihpc/multigrid.h
        scihpc/multigrid_base.cpp
        scihpc/multigrid_base.h
        )

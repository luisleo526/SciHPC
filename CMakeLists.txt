cmake_minimum_required(VERSION 3.21.0)
project(SciHPC)

if(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -O3")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -O3")
endif()

set(CMAKE_CXX_STANDARD 17)

file(GLOB SOURCES "scihpc/*.h" "scihpc/*.cpp")

add_executable(
    SciHPC
        ${SOURCES}
#        main.cpp
#        benchmarks_1d_wave.cpp
#        benchmarks_1d_IBE.cpp
#        benchmarks_2d_vortex.cpp
#        benchmarks_2d_zalesak.cpp
#        benchmarks_lsf_initialization.cpp
#        benchmarks_lsf_curvature_2d.cpp
#        benchmarks_lsf_curvature_3d.cpp
#        benchmarks_2d_cavity_flows.cpp
#        benchmarks_2d_dambreak.cpp
#        benchmarks_2d_RayleighTaylorInstability.cpp
        benchmarks_2d_bubble_rising.cpp
#        benchmarks_2d_multigrid_Poisson_Dirichlet.cpp
#        benchmarks_2d_multigrid_Poisson_Neumann.cpp
#        benchmarks_3d_dambreak.cpp
#        benchmarks_3d_multigrid_Poisson.cpp
)

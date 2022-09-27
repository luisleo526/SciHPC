cmake_minimum_required(VERSION 3.23)
project(SciHPC)

set(CMAKE_CXX_STANDARD 23)

add_executable(SciHPC main.cpp data_array.cpp data_array.h global.h structured_grid.cpp structured_grid.h)

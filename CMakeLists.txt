cmake_minimum_required(VERSION 3.23)
project(SciHPC)

set(CMAKE_CXX_STANDARD 23)

add_executable(SciHPC main.cpp Data.cpp Data.h global.h)

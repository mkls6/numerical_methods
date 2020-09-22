cmake_minimum_required(VERSION 3.17)
project(lu_decomposition)

set(CMAKE_CXX_STANDARD 17)

add_executable(lu_decomposition main.cpp src/matrix.cpp hpp/matrix.hpp src/linear_algebra.cpp hpp/linear_algebra.hpp)
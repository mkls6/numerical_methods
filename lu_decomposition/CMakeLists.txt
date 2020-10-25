cmake_minimum_required(VERSION 3.17)
project(lu_decomposition)

set(CMAKE_CXX_STANDARD 17)

add_executable(lu_decomposition main.cpp
        ../src/matrix.cpp
        ../include/matrix.hpp
        ../src/linear_algebra.cpp
        ../include/linear_algebra.hpp ../nl_system_solver/main.cpp)
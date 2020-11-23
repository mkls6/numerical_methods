cmake_minimum_required(VERSION 3.17)
project(nl_system_solver)

set(CMAKE_CXX_STANDARD 17)

add_executable(nl_system_solver main.cpp
        ../src/matrix.cpp
        ../include/matrix.hpp
        ../src/linear_algebra.cpp
        ../include/linear_algebra.hpp ../numerical_integration/main.cpp ../src/computational_math.cpp ../include/computational_math.hpp)

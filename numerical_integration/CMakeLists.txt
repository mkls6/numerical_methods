cmake_minimum_required(VERSION 3.17)
project(numerical_integration)

set(CMAKE_CXX_STANDARD 17)

add_executable(numerical_integration main.cpp
               ../include/computational_math.hpp
               ../src/computational_math.cpp
               main.cpp)

cmake_minimum_required(VERSION 3.5)
project(complete)

set(CMAKE_CXX_STANDARD 11)

add_executable(complete main.cpp point.hpp point.cpp foundation.hpp foundation.cpp)

cmake_minimum_required(VERSION 3.14)
project(cf)

set(CMAKE_CXX_STANDARD 17)

add_executable(cf main.cpp ContinuedFraction.cpp ContinuedFraction.h)
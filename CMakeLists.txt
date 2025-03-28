cmake_minimum_required(VERSION 3.30)
project(integrals)

set(CMAKE_CXX_STANDARD 23)

add_executable(integrals main.cpp
        J0.h
        integral_methods.h)

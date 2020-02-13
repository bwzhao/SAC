cmake_minimum_required(VERSION 3.14)
project(20190926)

set(CMAKE_CXX_STANDARD 11)
include_directories(/usr/local/include)
link_directories(/usr/local/lib)

set(CMAKE_CXX_FLAGS "-Wall -larmadillo -DARMA_DONT_USE_WRAPPER -framework Accelerate ${CMAKE_CXX_FLAGS}")
add_executable(20190926
        main.cpp
        Class_SpecFunc.cpp
        Class_SpecFunc.h
        Class_Calculate.cpp
        Class_Calculate.h
        Config.h
        Class_Spectral_Measurement.h)
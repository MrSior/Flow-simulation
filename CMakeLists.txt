cmake_minimum_required(VERSION 3.27)
project(fluid_project)

set(CMAKE_CXX_STANDARD 20)

add_executable(fluid_project
        fluid.cpp
        bits/stdc++.h)
add_executable(my-fluid my-fluid.cpp
        Fixed.h)
#target_compile_definitions(my-fluid PRIVATE -DTEST=FLOAT\(2\))
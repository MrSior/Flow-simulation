cmake_minimum_required(VERSION 3.27)
project(fluid_project)

set(CMAKE_CXX_STANDARD 20)

#set(CMAKE_CXX_FLAGS "-O3")

add_executable(FluidEngine Fixed.h
        customArray.h
        vectorField.h
        common.h
        FluidEngine.h
        main.cpp
        engineFactory.h
        inputParser.h)



target_compile_definitions(FluidEngine PRIVATE
        TYPES=FLOAT,DOUBLE,FIXED\(32,9\),FAST_FIXED\(40,5\),FAST_FIXED\(50,16\)
        SIZES=S\(10,250\),S\(36,84\)
)
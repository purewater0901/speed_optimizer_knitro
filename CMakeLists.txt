######################################################################
#
#  Copyright (c) 2016-2019 by Artelys.
#
######################################################################

# make a directory build in KTR, go there and run “cmake ..”.
# may need to specify generator in Windows:
# e.g.: for Visual Studio 2013 on x64, run cmake -G "Visual Studio 12 Win64" ..
# need separate build folders for Win32 and Win64
# The makefile will be created there. It uses KNITRODIR environment variable.

cmake_minimum_required(VERSION 2.8.6)

project(KnitroExamples)

find_package(OpenCV REQUIRED)
find_package(OpenMP)
if (OPENMP_FOUND)
  set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

include_directories(include)
include_directories($ENV{KNITRODIR}/include)

set(KTR_RELEASE 1200)
set(KNITRO_LIB  libknitro${KTR_RELEASE}.so gfortran)

if (APPLE)
  set(KNITRO_LIB "libknitro.dylib")
endif()

if(WIN32)
  set(KNITRO_LIB "knitro${KTR_RELEASE}.lib")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc")
endif()

link_directories($ENV{KNITRODIR}/lib)

add_executable(ExSpeedOptimizer src/ExSpeedOptimizer.cpp)
target_link_libraries(ExSpeedOptimizer ${KNITRO_LIB} ${CMAKE_DL_LIBS} ${OpenCV_LIBS})

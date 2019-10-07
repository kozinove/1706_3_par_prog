cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

project(samples)

include(cmake/function.cmake)

IF (WIN32)
  set(MPI_LIBRARIES "C:/Program Files/Microsoft HPC Pack 2008 SDK/Lib/i386/msmpi.lib" CACHE STRING "Path to MPI libraries")
  set(MPI_INCLUDE_PATH "C:/Program Files/Microsoft HPC Pack 2008 SDK/include" CACHE PATH "Path to MPI includes")
  message(STATUS ${MPI_INCLUDE_PATH})
  include_directories(${MPI_INCLUDE_PATH})

ELSE()
  find_package(MPI)
  include_directories(${MPI_INCLUDE_PATH}) 
  message(STATUS "Using MPI ${MPI_LIBRARIES}")
ENDIF()


add_subdirectory(function_lib)
add_subdirectory(main)

option(BTEST "build test?" OFF)

if(BTEST)
    add_subdirectory(gtest)
    add_subdirectory(test)
endif()


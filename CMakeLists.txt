cmake_minimum_required(VERSION 3.19.4)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS "-std=c++11 -O3")

project (parmetirs_ghost_v4 CXX)

find_package(MPI REQUIRED)
add_definitions(-DOMPI_SKIP_MPICXX)


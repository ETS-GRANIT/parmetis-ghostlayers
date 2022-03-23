cmake_minimum_required(VERSION 2.8.9)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS "--std=c++11 -O3 -fopenmp")

project (parmetirs_ghost_v4 C CXX)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
include_directories(.)

add_executable(main src/memory.c src/read_write.hpp src/read_write.cpp src/module_parmetis.hpp src/module_parmetis.cpp src/main.cpp)

target_link_libraries(main ${MPI_LIBRARIES} parmetis metis cgns)

if(MPI_COMPILE_FLAGS)
  set_target_properties(main PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(main PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

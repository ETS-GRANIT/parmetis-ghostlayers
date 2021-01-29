cmake_minimum_required(VERSION 2.8.9)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS "-std=c++11 -O3")

project (parmetirs_ghost_v4 CXX)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

FIND_PATH(PARMETIS_INCLUDE_DIR parmetis.h
  /usr/local/include
  /usr/include
  )

FIND_LIBRARY(PARMETIS_LIBRARY parmetis
  /usr/local/lib
  /usr/lib
  )

FIND_LIBRARY(METIS_LIBRARY metis
  /usr/local/lib
  /usr/lib
  )

IF(PARMETIS_INCLUDE_DIR)
  IF(PARMETIS_LIBRARY)
    SET( PARMETIS_LIBRARIES ${PARMETIS_LIBRARY} ${METIS_LIBRARY})
    SET( PARMETIS_FOUND "YES" )
  ENDIF(PARMETIS_LIBRARY)
ENDIF(PARMETIS_INCLUDE_DIR)

file(GLOB SOURCES "*.cpp")
add_executable(main ${SOURCES})

target_link_libraries(main ${MPI_LIBRARIES} ${PARMETIS_LIBRARY} ${METIS_LIBRARY})

if(MPI_COMPILE_FLAGS)
  set_target_properties(main PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(main PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

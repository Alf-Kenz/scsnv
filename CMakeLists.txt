cmake_minimum_required (VERSION 3.1)
project (scsnv C CXX)

if (POLICY CMP0025)
  cmake_policy(SET CMP0025 NEW)
endif ()

set (CMAKE_CXX_STANDARD 11)
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )

set(DEBUG_FLAGS -g)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DEBUG_FLAGS} -Wall -pedantic -Wextra -fno-omit-frame-pointer -g")
#set(HDF5_USE_STATIC_LIBRARIES ON)

find_package( HDF5 REQUIRED COMPONENTS C CXX HL)
INCLUDE_DIRECTORIES(SYSTEM ${HDF5_INCLUDE_DIRS})
MESSAGE( "HDF5 LIBRARY: ${HDF5_LIBRARIES}" )
MESSAGE( "HDF5 LIBRARY CXX: ${HDF5_CXX_LIBRARIES}" )
MESSAGE( "HDF5 INCLUDE: ${HDF5_INCLUDE_DIRS}" )

MESSAGE( "BUILD TYPE: ${CMAKE_BUILD_TYPE}" )
MESSAGE( "INSTALL PREFIX: ${CMAKE_INSTALL_PREFIX}" )

SET( EXTERNAL_DIR "${PROJECT_SOURCE_DIR}/external" )

#add_subdirectory("bwa")
add_subdirectory("src")

# mamba install -c conda-forge clang clangxx

# export CXX=/home/jlanglie/anaconda3/envs/scsnvpy/bin/clang++; export CC=/home/jlanglie/anaconda3/envs/scsnvpy/bin/clang
# export CXX=; export CC=
# cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-fsanitize=address  " ..

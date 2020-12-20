project(HODLRlib)
cmake_minimum_required(VERSION 3.0)
set(CMAKE_CXX_STANDARD 11)

################ FINDING MKL ####################

# Checking for MKL:
# Currently disabling MKL. Sometimes a linker error pops up
# TODO: Look into this
#if (DEFINED ENV{MKLROOT})
#    add_compile_definitions(MKL_ENABLED)
#    set (MKLROOT "$ENV{MKLROOT}")
#    set (LDFLAGS "${LDFLAGS} -DMKL_ENABLED -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group")
#    set (MKL_FLAGS "-DMKL_LP64 -DMKL_ENABLED -I${MKLROOT}/include")
#endif ()

################ FINDING EIGEN ####################

# Change this accordingly to point to the Eigen root folder
# By default, it will pick up the environment variable EIGEN_PATH
set(EIGEN_PATH "$ENV{EIGEN_PATH}")

if(NOT EIGEN_PATH)
    message ( FATAL_ERROR "Please point the environment variable EIGEN_PATH to the root directory of your Eigen3 installation.")
endif()

################ SETTING DATATYPE TO USE:FLOAT, DOUBLE AND COMPLEX ####################

# Set this to: "float", "double", "complex32" and "complex64"
set(DTYPE
    "double" CACHE STRING
    "Matrix entry type: \"float\", \"double\", \"complex32\" or \"complex64\""
   )
set_property(CACHE DTYPE PROPERTY STRINGS float double complex32 complex64)

if(DTYPE STREQUAL "float")
    set(DTYPE_FLAG "-DUSE_FLOAT")
elseif(DTYPE STREQUAL "double")
    set(DTYPE_FLAG "-DUSE_DOUBLE")
elseif(DTYPE STREQUAL "complex32")
    set(DTYPE_FLAG "-DUSE_COMPLEX32")
elseif(DTYPE STREQUAL "complex64")
    set(DTYPE_FLAG "-DUSE_COMPLEX64")
else()
    message(FATAL_ERROR "Invalid Choice for Datatype. Exiting...")
endif()

################ COMPILERS: APPROPRIATELY USE INTEL OR GNU ####################

# Toggle this to 1 when Intel Compilers are available:
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DTYPE_FLAG} ${MKL_FLAGS} --std=c++11 -O3 -ipo -qopenmp -xHost -funroll-loops -ffast-math -ffinite-math-only -g -liomp5 -lpthread -ldl")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DTYPE_FLAG} ${MKL_FLAGS} --std=c++11 -O3 -funroll-loops -ffast-math -ffinite-math-only -g -fopenmp -lpthread -ldl")
endif ()

# Flags for testing coverage information
# Need by coveralls:
set(CMAKE_CXX_FLAGS_COVERAGE "${CMAKE_CXX_FLAGS} --coverage")
set(CMAKE_EXE_LINKER_FLAGS_COVERAGE "${CMAKE_EXE_LINKER_FLAGS} --coverage")
set(CMAKE_SHARED_LINKER_FLAGS_COVERAGE "${CMAKE_SHARED_LINKER_FLAGS} --coverage")

################ ADDING OF THE NECESSARY HEADER FOLDERS ####################
# Adding folders which contain necessary header files:
include_directories("${EIGEN_PATH}"
                    "include"
                   )

################ MENTIONING C++ FILES TO LINK WITH HEADERS + FILE TO BUILD ####################

set(HODLRlib_SRCS
    ${PROJECT_SOURCE_DIR}/src/HODLR_Matrix.cpp
    ${PROJECT_SOURCE_DIR}/src/LowRank.cpp
    ${PROJECT_SOURCE_DIR}/src/HODLR_Node.cpp
    ${PROJECT_SOURCE_DIR}/src/HODLR.cpp
    ${PROJECT_SOURCE_DIR}/src/HODLR_NonSPD.cpp
    ${PROJECT_SOURCE_DIR}/src/HODLR_SPD.cpp
    ${PROJECT_SOURCE_DIR}/src/KDTree.cpp
   )

add_library(HODLR ${HODLRlib_SRCS})

option(HODLR_BUILD_TESTS "Build tests for the HODLR library" ON)
if(NOT("${CMAKE_BUILD_TYPE}" STREQUAL "COVERAGE") AND ${HODLR_BUILD_TESTS})
    add_subdirectory(examples)
    add_subdirectory(test)
endif()

################ SUMMARY OF BUILD INFORMATION ################

message("")
message("#################")
message("# Build options #")
message("#################")
message("")
message("Input File       : " ${INPUT_FILE})
message("Data Type        : " ${DTYPE})
message("CXX Compiler     : " ${CMAKE_CXX_COMPILER})
message("CXX_FLAGS        : " ${CMAKE_CXX_FLAGS})
message("Eigen dir        : " ${EIGEN_PATH})
# message("MKL dir          : " ${MKLROOT})
message("Build dir        : " ${CMAKE_BINARY_DIR})
message("")

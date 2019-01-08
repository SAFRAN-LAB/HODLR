project (HODLRlib)
cmake_minimum_required (VERSION 3.0)
set (CMAKE_CXX_STANDARD 11)
# Provide the input file that you want to solve for:
set (INPUT_FILE "./test_HODLR.cpp")
# Here set the name for the output executable:
set (OUTPUT_EXECUTABLE "test_HODLR")
# Set this to: "float", "double", "complex32" and "complex64"
set (DTYPE "double")
# Change this accordingly to point to the Eigen root folder
# By default, it will pick up the environment variable EIGEN_PATH
set(EIGEN_PATH "$ENV{EIGEN_PATH}")

################ SETTING DATATYPE TO USE:FLOAT, DOUBLE AND COMPLEX ####################

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

############### SETTING COMPILER FLAGS ##################################

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DTYPE_FLAG} -O3 -funroll-loops -ffast-math -ffinite-math-only -g --coverage  -fopenmp -lpthread -ldl")

################ ADDING OF THE NECESSARY HEADER FOLDERS ####################

# Adding folders which contain necessary header files:
include_directories("${EIGEN_PATH}" 
                    "../src"
                   )

################ MENTIONING C++ FILES TO LINK WITH HEADERS + FILE TO BUILD ####################
add_executable(${OUTPUT_EXECUTABLE} 
               ${INPUT_FILE}
               ../src/HODLR_Matrix.cpp
               ../src/HODLR_Node.cpp
               ../src/HODLR_Tree.cpp
               ../src/HODLR_Tree_NonSPD.cpp
               ../src/HODLR_Tree_SPD.cpp
               ../src/KDTree.cpp
              )

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
message("Build dir        : " ${CMAKE_BINARY_DIR})
message("")

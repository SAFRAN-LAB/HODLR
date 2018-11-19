project(tryingOut)
cmake_minimum_required (VERSION 3.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -fopenmp")

# Change this accordingly to point to the Eigen root folder:
SET(EIGEN_PATH "$ENV{EIGEN_PATH}")
IF(NOT EIGEN_PATH)
    MESSAGE( FATAL_ERROR "Please point the environment variable EIGEN_PATH to the root directory of your Eigen3 installation.")
ENDIF()

# Adding folders which contain necessary header files:
INCLUDE_DIRECTORIES("${EIGEN_PATH}" "src")

add_executable(HODLR 
               examples/testHODLR.cpp
               src/HODLR_Matrix.cpp
               src/HODLR_Node.cpp
               src/HODLR_Tree.cpp
              )

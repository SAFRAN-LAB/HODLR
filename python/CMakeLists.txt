cmake_minimum_required (VERSION 3.0)
project(HODLRlib_python)
add_subdirectory(pybind11)

# Change this accordingly to point to the Eigen root folder
# By default, it will pick up the environment variable EIGEN_PATH
set(EIGEN_PATH "$ENV{EIGEN_PATH}")
set(HODLR_PATH "$ENV{HODLR_PATH}")

set(DTYPE_FLAG "-DUSE_DOUBLE")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DTYPE_FLAG} --std=c++11 -O3 -funroll-loops -ffast-math -ffinite-math-only")

################ ADDING OF THE NECESSARY HEADER FOLDERS ####################
# Adding folders which contain necessary header files:
include_directories("${EIGEN_PATH}" 
                    "${HODLR_PATH}/include"
                   )

################ MENTIONING C++ FILES TO LINK WITH HEADERS + FILE TO BUILD ####################
set(HODLRlib_SRCS
    ${HODLR_PATH}/src/HODLR_Matrix.cpp
    ${HODLR_PATH}/src/KDTree.cpp
    ${HODLR_PATH}/src/Matrix_Factorizer.cpp
    ${HODLR_PATH}/src/HODLR_Node.cpp
    ${HODLR_PATH}/src/HODLR_Tree.cpp
    ${HODLR_PATH}/src/HODLR_Tree_NonSPD.cpp
    ${HODLR_PATH}/src/HODLR_Tree_SPD.cpp
   )

pybind11_add_module(pyhodlrlib
                    pyhodlrlib.cpp
                    ${HODLRlib_SRCS}
                   )

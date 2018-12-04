project (HODLR)
cmake_minimum_required (VERSION 3.0)
set (CMAKE_CXX_STANDARD 11)

# Checking for MKL:
if (DEFINED ENV{MKLROOT})
    message("Found MKL!")
    set (MKLROOT "$ENV{MKLROOT}")
    set (LDFLAGS "-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group")
    set (MKL_FLAGS "-DMKL_LP64 -I${MKLROOT}/include")
endif ()

# Toggle this to 1 when Intel Compilers are available:
if (0)
    message("Using Intel Compilers!")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MKL_FLAGS} -O3 -xHost -ipo -fopenmp -ffast-math -g -liomp5 -lpthread -ldl")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MKL_FLAGS} -O3 -ffast-math -g -fopenmp -lpthread -ldl")
endif ()

# Change this accordingly to point to the Eigen root folder:
set(EIGEN_PATH "$ENV{EIGEN_PATH}")
if(NOT EIGEN_PATH)
    message ( FATAL_ERROR "Please point the environment variable EIGEN_PATH to the root directory of your Eigen3 installation.")
endif()

# Adding folders which contain necessary header files:
include_directories("${EIGEN_PATH}" "src")

add_executable(example 
               examples/testHODLR.cpp
               src/HODLR_Matrix.cpp
               src/HODLR_Node.cpp
               src/HODLR_Tree.cpp
               src/KDTree.cpp
              )

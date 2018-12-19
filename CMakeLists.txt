project (HODLR)
cmake_minimum_required (VERSION 3.0)
set (CMAKE_CXX_STANDARD 11)

################ FINDING MKL ####################

# Checking for MKL:
if (DEFINED ENV{MKLROOT})
    add_compile_definitions(MKL_ENABLED)
    set (MKLROOT "$ENV{MKLROOT}")
    set (LDFLAGS "${LDFLAGS} -DMKL_ENABLED -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group")
    set (MKL_FLAGS "-DMKL_LP64 -DMKL_ENABLED -I${MKLROOT}/include")
endif ()

################ FINDING EIGEN ####################

# Change this accordingly to point to the Eigen root folder:
set(EIGEN_PATH "$ENV{EIGEN_PATH}")
if(NOT EIGEN_PATH)
    message ( FATAL_ERROR "Please point the environment variable EIGEN_PATH to the root directory of your Eigen3 installation.")
endif()

################ COMPILERS: APPROPRIATELY USE INTEL OR GNU ####################

# Toggle this to 1 when Intel Compilers are available:
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MKL_FLAGS} -O3 -ipo -qopenmp -xHost -funroll-loops -ffast-math -ffinite-math-only -g -liomp5 -lpthread -ldl")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MKL_FLAGS} -O3 -funroll-loops -ffast-math -ffinite-math-only -g -fopenmp -lpthread -ldl")
endif ()

################ ADDING OF THE NECESSARY HEADER FOLDERS ####################
# Adding folders which contain necessary header files:
include_directories("${EIGEN_PATH}" 
                    "src"
                   )

################ MENTION C++ FILES TO LINK WITH HEADERS + FILE TO BUILD ####################
add_executable(example 
               examples/testHODLR.cpp
               src/HODLR_Matrix.cpp
               src/HODLR_Node.cpp
               src/HODLR_Tree.cpp
               src/KDTree.cpp
              )

################ SUMMARY OF BUILD INFORMATION ################

message("")
message("#################")
message("# Build options #")
message("#################")
message("")
message("CXX Compiler     : " ${CMAKE_CXX_COMPILER})
message("CXX_FLAGS        : " ${CMAKE_CXX_FLAGS})
message("Eigen dir        : " ${EIGEN_PATH})
message("MKL dir          : " ${MKLROOT})
message("Build dir        : " ${CMAKE_BINARY_DIR})
message("")

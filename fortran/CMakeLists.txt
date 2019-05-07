project (HODLRlib_fortran)
cmake_minimum_required (VERSION 3.0)
set (CMAKE_CXX_STANDARD 11)
enable_language(Fortran)

# Provide the input file that you want to solve for:
# This can also be passed as a command line argument:
if(NOT INPUT_FILE)
    set (INPUT_FILE "example/example.f90")
endif()
# Here set the name for the output executable:
# This can also be passed as a command line argument:
if(NOT OUTPUT_EXECUTABLE)
    set (OUTPUT_EXECUTABLE "example")
endif()

if("${CMAKE_BUILD_TYPE}" STREQUAL "TESTING")
    set (INPUT_FILE "test/test_HODLR.f90")
    set (OUTPUT_EXECUTABLE "test_HODLR")
endif()

set(EIGEN_PATH "$ENV{EIGEN_PATH}")
set(HODLR_PATH "$ENV{HODLR_PATH}")

set(DTYPE_FLAG "-DUSE_DOUBLE")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DTYPE_FLAG} --std=c++11 -O3 -funroll-loops -ffast-math -ffinite-math-only")
set(CMAKE_Fortran_FLAGS "-lstdc++")

################ ADDING OF THE NECESSARY HEADER FOLDERS ####################
# Adding folders which contain necessary header files:
include_directories("${EIGEN_PATH}" 
                    "${HODLR_PATH}/include"
                    "${HODLR_PATH}/fortran/wrapper"
                    "${HODLR_PATH}/fortran/example"
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

add_executable(${OUTPUT_EXECUTABLE}
               ${INPUT_FILE}
               ${HODLRlib_SRCS}
               "${HODLR_PATH}/fortran/wrapper/fortran_wrapper.cpp"
               "${HODLR_PATH}/fortran/wrapper/hodlr_mod.f90"
              )

set_property(TARGET ${OUTPUT_EXECUTABLE} PROPERTY LINKER_LANGUAGE Fortran)

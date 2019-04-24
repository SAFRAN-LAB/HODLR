project (HODLRlib_fortran)
cmake_minimum_required (VERSION 3.0)
set (CMAKE_CXX_STANDARD 11)
enable_language(Fortran)

# Provide the input file that you want to solve for:
# This can also be passed as a command line argument:
if(NOT INPUT_FILE)
    set (INPUT_FILE "example/test.f90")
endif()
# Here set the name for the output executable:
# This can also be passed as a command line argument:
if(NOT OUTPUT_EXECUTABLE)
    set (OUTPUT_EXECUTABLE "test")
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
                    "${HODLR_PATH}/fortran/kernel"
                    "${HODLR_PATH}/fortran/example"
                   )

################ MENTIONING C++ FILES TO LINK WITH HEADERS + FILE TO BUILD ####################

set(HODLRlib_SRCS
    ${HODLR_PATH}/src/HODLR_Matrix.cpp
   )

set(kernellib_SRCS
    ${HODLR_PATH}/fortran/kernel/kernel_capi.cpp
    ${HODLR_PATH}/fortran/kernel/kernel_mod.f90
   )

add_executable(${OUTPUT_EXECUTABLE}
               ${INPUT_FILE}
               ${HODLRlib_SRCS}
               ${kernellib_SRCS}
               "${HODLR_PATH}/fortran/example/Kernel.cpp"
              )

set_property(TARGET test PROPERTY LINKER_LANGUAGE Fortran)

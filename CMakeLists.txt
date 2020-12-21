cmake_policy(SET CMP0048 NEW)
project(HODLRlib VERSION 3.1415)
cmake_minimum_required(VERSION 3.0)
set(CMAKE_CXX_STANDARD 11)

################ OPTIONS ####################

option(HODLR_BUILD_TESTS "Build tests for the HODLR library" ON)
option(USE_MKL "Use Intel MKL with Eigen." OFF)

################ FIND DEPENDENCIES ####################

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

find_package(OpenMP REQUIRED)

if(${USE_MKL})
  add_compile_definitions(MKL_ENABLED)
  set(BLA_VENDOR Intel10_64lp)
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
endif ()

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(find_or_download)

find_or_download(Eigen3)

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
set(CMAKE_CXX_FLAGS)
add_compile_options(
  --std=c++11
  $<$<CONFIG:>:-O3>
  -funroll-loops -ffast-math -ffinite-math-only
)
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ipo -xHost -liomp5")
endif ()

# Flags for testing coverage information
# Need by coveralls:
set(CMAKE_CXX_FLAGS_COVERAGE "${CMAKE_CXX_FLAGS} --coverage")
set(CMAKE_EXE_LINKER_FLAGS_COVERAGE "${CMAKE_EXE_LINKER_FLAGS} --coverage")
set(CMAKE_SHARED_LINKER_FLAGS_COVERAGE "${CMAKE_SHARED_LINKER_FLAGS} --coverage")

################ MAIN BUILD INSTRUCTIONS ####################

add_library(HODLR
    ${PROJECT_SOURCE_DIR}/src/HODLR_Matrix.cpp
    ${PROJECT_SOURCE_DIR}/src/LowRank.cpp
    ${PROJECT_SOURCE_DIR}/src/HODLR_Node.cpp
    ${PROJECT_SOURCE_DIR}/src/HODLR.cpp
    ${PROJECT_SOURCE_DIR}/src/HODLR_NonSPD.cpp
    ${PROJECT_SOURCE_DIR}/src/HODLR_SPD.cpp
    ${PROJECT_SOURCE_DIR}/src/KDTree.cpp
   )
target_link_libraries(HODLR
  PUBLIC Eigen3::Eigen $<$<BOOL:USE_MKL>:${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}>
  PRIVATE Threads::Threads OpenMP::OpenMP_CXX
)
target_include_directories(HODLR
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_compile_definitions(HODLR PUBLIC ${DTYPE_FLAG})

if(${HODLR_BUILD_TESTS})
    add_subdirectory(examples)
    add_subdirectory(test)
endif()

################ SUMMARY OF BUILD INFORMATION ################

message("")
message("#################")
message("# Build options #")
message("#################")
message("")
message("Data Type        : " ${DTYPE})
message("CXX Compiler     : " ${CMAKE_CXX_COMPILER})
message("CXX_FLAGS        : " ${CMAKE_CXX_FLAGS})
message("Eigen dir        : " ${EIGEN3_ROOT_DIR})
message("MKL dir          : " $ENV{MKLROOT})
message("Build dir        : " ${CMAKE_BINARY_DIR})
message("")

################ Install instruction ################
# Create version file for cmake package
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  HODLRConfigVersion.cmake
  VERSION ${PACKAGE_VERSION}
  COMPATIBILITY SameMajorVersion
)
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/HODLRConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/cmake/HODLR
)
# Create config file of HODLR
# No public dependencies so a Config file with the targets is enough
install(TARGETS HODLR EXPORT HODLRTargets DESTINATION lib)
install(EXPORT HODLRTargets
  FILE HODLRConfig.cmake
  NAMESPACE HODLR::
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/cmake/HODLR
)
# Copy headers to install directory
install(
  DIRECTORY ${CMAKE_SOURCE_DIR}/include/HODLR
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include
)

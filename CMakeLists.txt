project (HODLR C CXX)
cmake_minimum_required (VERSION 2.8)
set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# Set the version.
set (HODLR_VERSION_MAJOR 3)
set (HODLR_VERSION_MINOR 14)
set (HODLR_VERSION_PATCH 0)
set (HODLR_VERSION
    ${HODLR_VERSION_MAJOR}.${HODLR_VERSION_MINOR}.${HODLR_VERSION_PATCH})
set (HODLR_ABI_VERSION 3.14.0)

# Configure the path structure.
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set (CMAKE_RELEASE_POSTFIX "")
set (CMAKE_DEBUG_POSTFIX "-debug")

# Turn on testing.
ENABLE_TESTING()

# Find the dependencies.
find_package (Eigen REQUIRED)
if (EIGEN_FOUND)
    message ("-- Found Eigen version ${EIGEN_VERSION}: ${EIGEN_INCLUDE_DIRS}")
endif (EIGEN_FOUND)

# Add the include directories.
include_directories (header ${EIGEN_INCLUDE_DIRS})

# Install the public headers.
file (GLOB HODLR_HDRS ${CMAKE_SOURCE_DIR}/header/*.hpp)
install (FILES ${HODLR_HDRS} DESTINATION include)

# Heuristic for determining LIB_SUFFIX. FHS recommends that 64-bit systems
# install native libraries to lib64 rather than lib. Most distros seem to
# follow this convention with a couple notable exceptions (Debian-based and
# Arch-based distros) which we try to detect here.
if (CMAKE_SYSTEM_NAME MATCHES "Linux" AND
        NOT DEFINED LIB_SUFFIX AND
        NOT CMAKE_CROSSCOMPILING AND
        CMAKE_SIZEOF_VOID_P EQUAL "8" AND
        NOT EXISTS "/etc/debian_version" AND
        NOT EXISTS "/etc/arch-release")
    set (LIB_SUFFIX "64")
endif ()

# Set up the compiler flags and optimizations.
set (HODLR_FLAGS "${CMAKE_CXX_FLAGS} -DNDEBUG -ffast-math")

if (CMAKE_COMPILER_IS_GNUCC)
    # Linux
    if (CMAKE_SYSTEM_NAME MATCHES "Linux")
        if (NOT GCC_VERSION VERSION_LESS 4.2)
            set (HODLR_FLAGS "${HODLR_FLAGS} -march=native -mtune=native")
        endif (NOT GCC_VERSION VERSION_LESS 4.2)
    endif (CMAKE_SYSTEM_NAME MATCHES "Linux")
    # Mac OS X
    if (CMAKE_SYSTEM_NAME MATCHES "Darwin")
        set (HODLR_FLAGS "${HODLR_FLAGS} -msse3")
        # Use of -fast only applicable for Apple's GCC
        # Assume this is being used if GCC version < 4.3 on OSX
        execute_process (COMMAND ${CMAKE_C_COMPILER}
                         ARGS ${CMAKE_C_COMPILER_ARG1} -dumpversion
                         OUTPUT_VARIABLE GCC_VERSION
                         OUTPUT_STRIP_TRAILING_WHITESPACE)
        if (GCC_VERSION VERSION_LESS 4.3)
            set (HODLR_FLAGS "${HODLR_FLAGS} -fast")
        endif (GCC_VERSION VERSION_LESS 4.3)
    endif (CMAKE_SYSTEM_NAME MATCHES "Darwin")
endif (CMAKE_COMPILER_IS_GNUCC)

if (CMAKE_C_COMPILER_ID STREQUAL "Clang")
    # Use of -O4 requires use of gold linker & LLVM-gold plugin, which might
    # well not be present / in use and without which files will compile, but
    # not link ('file not recognized') so explicitly check for support
    include (CheckCCompilerFlag)
    check_c_compiler_flag ("-O4" HAVE_LTO_SUPPORT)
    if (HAVE_LTO_SUPPORT)
        message (STATUS "Enabling link-time optimization (-O4)")
        set (HODLR_FLAGS "${HODLR_FLAGS} -O4")
    else ()
        message (STATUS
                 "Compiler/linker does not support link-time optimization (-O4), disabling.")
    endif (HAVE_LTO_SUPPORT)
endif ()

# GCC is not strict enough by default, so enable most of the warnings.
if (UNIX)
    set (HODLR_FLAGS
        "${HODLR_FLAGS} -Werror -Wall -Wextra -Wno-unknown-pragmas -Wno-sign-compare -Wno-unused-parameter -Wno-missing-field-initializers")
endif (UNIX)

# Compile the library.
file (GLOB HODLR_SRC ${CMAKE_SOURCE_DIR}/src/*.cpp)
set (CMAKE_CXX_FLAGS "${HODLR_FLAGS}")
add_library (hodlr STATIC ${HODLR_SRC})
install (TARGETS hodlr
         RUNTIME DESTINATION bin
         LIBRARY DESTINATION lib${LIB_SUFFIX}
         ARCHIVE DESTINATION lib${LIB_SUFFIX})

# Build the Kepler demo.
add_executable (kepler examples/kepler_test.cpp)
target_link_libraries (kepler hodlr)

# Build the tests.
macro (HODLR_TEST NAME KERNEL DIM)
    add_executable (${NAME}_test examples/HODLR_Test.cpp)
    target_link_libraries (${NAME}_test hodlr)
    set (CMAKE_CXX_FLAGS "${HODLR_FLAGS} -D${KERNEL} -D${DIM}")
    add_test (NAME ${NAME}_test
              COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${NAME}_test)
endmacro (HODLR_TEST)

hodlr_test (oned_gaussian GAUSSIAN ONE)
hodlr_test (twod_gaussian GAUSSIAN TWO)

macro(find_or_download PACKAGE)
  set(DEPENDENCY_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/dependencies)
  # Update search path and use regular find_package to add dependency
  find_package(
    ${PACKAGE} QUIET
    HINTS ${CMAKE_SOURCE_DIR}/dependencies ${CMAKE_INSTALL_PREFIX}
  )

  if(${${PACKAGE}_FOUND})
    message(STATUS "Found dependency ${PACKAGE} installed in system.")
  else()
    message(STATUS "Suitable version of ${PACKAGE} not found in system.")
    message(STATUS "Downloading ${PACKAGE} and building from source.")
    # Prepare download instructions for dependency
    configure_file(
      ${CMAKE_SOURCE_DIR}/cmake/${PACKAGE}_download.cmake.in
      ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE}_download/CMakeLists.txt
      @ONLY
    )

    # Configure step for download instructions
    execute_process(
      COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
      RESULT_VARIABLE result
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE}_download
      OUTPUT_QUIET
    )
    if(result)
      message(FATAL_ERROR "Download of dependency ${PACKAGE} failed: ${result}")
    endif()

    # Download, build and install dependency according to instructions
    execute_process(
      COMMAND ${CMAKE_COMMAND} --build .
      RESULT_VARIABLE result
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE}_download
    )
    if(result)
      message(FATAL_ERROR "Build of dependency ${PACKAGE} failed: ${result}.")
    endif()
    # Update search path and use regular find_package to add dependency
    find_package(${PACKAGE}
      REQUIRED NO_DEFAULT_PATH
      PATHS "${DEPENDENCY_INSTALL_PREFIX}"
    )
    message(STATUS "Using ${PACKAGE} from ${DEPENDENCY_INSTALL_PREFIX}.")
  endif()
endmacro()

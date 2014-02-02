find_program (PYTHON "python")

if (PYTHON)
    message ("-- Found Python interpreter: ${PYTHON}")

    set (SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
    set (SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
    # set (INIT_PY_IN  "${CMAKE_CURRENT_SOURCE_DIR}/george/__init__.py.in")
    # set (INIT_PY     "${CMAKE_CURRENT_SOURCE_DIR}/george/__init__.py")
    set (DEPS        "${CMAKE_CURRENT_SOURCE_DIR}/hodlr/*")
    set (OUTPUT      "${CMAKE_CURRENT_BINARY_DIR}/build/timestamp")

    configure_file (${SETUP_PY_IN} ${SETUP_PY})
    # configure_file (${INIT_PY_IN} ${INIT_PY})

    add_custom_command (OUTPUT ${OUTPUT}
                        COMMAND ${PYTHON} ${SETUP_PY} build
                        COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
                        DEPENDS ${DEPS})

    add_custom_target (python ALL DEPENDS ${OUTPUT})

    install (CODE
             "EXECUTE_PROCESS(COMMAND ${PYTHON} ${SETUP_PY} install ${PYTHON_INSTALL_OPTIONS})")

else (PYTHON)
    message ("-- Did not find Python interpreter")
endif ()

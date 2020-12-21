foreach(TARGET
        tutorial
        example_lowrank_factorization
        example
        example_RPY
        example_matern
        example_conducting_sphere_plate
       )

  add_executable(${TARGET} ${TARGET}.cpp)
  target_link_libraries(${TARGET} HODLR)
endforeach()

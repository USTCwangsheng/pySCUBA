foreach(program testtransform testcrdtree testcalculators)
   add_executable(${program} ${program}.cpp)
   target_link_libraries(${program} ${ALL_LIBS})
endforeach()

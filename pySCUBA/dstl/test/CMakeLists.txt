foreach(program testtree testnbstree testpermutatioin testparetofront testnestedcontainers testtopn testsortidx)
   add_executable(${program} ${program}.cpp)
   target_link_libraries(${program} ${ALL_LIBS})
endforeach()

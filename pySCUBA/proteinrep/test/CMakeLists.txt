foreach(program testresiduestate testpocketmodeler testisolatedpocket getsequencefrompdb)
   add_executable(${program} ${program}.cpp)
   target_link_libraries(${program} ${ALL_LIBS})
endforeach()

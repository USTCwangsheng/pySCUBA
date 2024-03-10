set(PROGRAMS findssas extractligand testsmartsfragment testatomtypes analyze_iter contactstatistics calcrefdistr getcontactligandnames)
foreach(program ${PROGRAMS})
 add_executable(${program} ${program}.cpp)
 target_link_libraries(${program} combine openbabel cliquer)
endforeach(program)


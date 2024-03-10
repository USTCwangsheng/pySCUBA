foreach(program hbgeometry testbackbonesasa loophbsasa testbackboneloop nativevsdecoyloops testmainchain checksitesfile gentrialloops testcoresegments optimizelinkedchain decomposesiteenergies testenergyfunctions extractlooptochain genloop testalign testscorepep testbackboneenergy testminimizer gaminimizer testbackbonebuilder testchainpackminimizer backbonedecoys testtvscore)   add_executable(${program} ${program}.cpp)
   target_link_libraries(${program}  ${ALL_LIBS})
endforeach()

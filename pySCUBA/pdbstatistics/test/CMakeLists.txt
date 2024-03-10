foreach(program estimatephipsidistr testphipsidistr extractphipsi linkerlengths looplocations testscoresitepair estimatepbtetrabase testproteinblock testssclassifier testvirtualsidechain estimatecbpairene)
   add_executable(${program} ${program}.cpp)
   target_link_libraries(${program} ${ALL_LIBS})
endforeach()

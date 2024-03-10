foreach(program testintrctmol mkparfiletemplates testmolmodeler genrandomconfig randomchain SCUBAChangeResidues buildsidechains preprocesspdb SCUBASketch comparess extendloopread reduceloopread readss makeLVGintendedss)
   add_executable(${program} ${program}.cpp)
   target_link_libraries(${program} ${ALL_LIBS})
endforeach()

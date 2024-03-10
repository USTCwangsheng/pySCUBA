foreach(program printpartemp SCUBASD testloopsampler readsampledloops testgroupdistance testligandpep testligpsampler readligpconfig SCUBALoopSampler SCUBALoopReader SCUBAHBScreen)
   add_executable(${program} ${program}.cpp)
   target_link_libraries(${program} ${ALL_LIBS})
endforeach()



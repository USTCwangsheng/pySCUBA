link_directories(${PROJECT_BINARY_DIR}/lib)



add_executable(DesignSeq TestDesignSeq.cpp)

target_link_libraries (DesignSeq ${ALL_LIBS})


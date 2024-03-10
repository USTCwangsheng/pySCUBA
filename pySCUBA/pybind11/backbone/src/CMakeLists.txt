set(BACKBONE_SOURCE_FILES pb_backbone.cpp)
add_library(backbone ${BACKBONE_SOURCE_FILES})
#target_link_libraries(backbone ${PythonLibs_LIBRARIES} ${Boost_LIBRARIES})
target_link_libraries(backbone ${PythonLibs_LIBRARIES})
set_target_properties(backbone PROPERTIES PREFIX "")
#file(GLOB USER_CPP ${SCUBAPATH}/*/src/*.cpp)

#set(PRO_SOURCE_FILES ${SCUBAPATH}/backbone/src/backbonesite.cpp ${SCUBAPATH}/backbone/src/backbonebuilder.cpp ${SCUBAPATH}/backbone/src/closealoop.cpp ${SCUBAPATH}/geometry/src/calculators.cpp ${SCUBAPATH}/geometry/src/localframe.cpp ${SCUBAPATH}/proteinrep/src/idealgeometries.cpp ${SCUBAPATH}/proteinrep/src/pdbreader.cpp ${SCUBAPATH}/pdbstatistics/src/phipsidistr.cpp ${SCUBAPATH}/dataio/src/datapaths.cpp ${SCUBAPATH}/geometry/src/quaternioncrd.cpp ${SCUBAPATH}/dataio/src/inputlines.cpp ${SCUBAPATH}/dataio/src/splitstring.cpp ${SCUBAPATH}/loopclosure/src/loopsolver.cpp ${SCUBAPATH}/loopclosure/src/tripepclosure.cpp ${SCUBAPATH}/geometry/src/rigidtransform.cpp ${SCUBAPATH}/geometry/src/quatfit.cpp)

pybind11_add_module(pb_backbone pb_backbone.cpp)
target_link_libraries(pb_backbone PRIVATE ${ALL_LIBS})
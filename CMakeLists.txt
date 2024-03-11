cmake_minimum_required(VERSION 3.7)

project(pySCUBA)
set(MY_ROOT ${PROJECT_SOURCE_DIR})
message(${MY_ROOT})
# Find necessary packages
find_package(PythonLibs REQUIRED)
set(PythonLibs_LIBRARIES ~/anaconda3/envs/pySCUBA/lib/libpython3.8.so)
message("python libraries: ${PythonLibs_LIBRARIES}")
find_package(LAPACK REQUIRED)
set(LAPACK_LIBRARIES ~/.local/lapack/liblapacke.so)
message("lapack libraries: ${LAPACK_LIBRARIES}")
include_directories(${PYTHON_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/pySCUBA ${CMAKE_SOURCE_DIR}/pySCUBA/include/)
message("PythonInclude:${PYTHON_INCLUDE_DIR}")



add_definitions(-DLIBARCHIVE_STATIC -fPIC)
set(pybind11_INCLUDE_DIRS ./3rdparty/pybind11)
include_directories(./3rdparty/pybind11/include)

set(CMAKE_CXX_STANDARD 11)
add_compile_options("-Wno-return-type")
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DDEBUG_LEVEL=0 -O2 ")

set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/pySCUBA/pybind11/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/pySCUBA/pybind11/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/pySCUBA/cpp_bin)
set(SCUBA_LIBS nspsampling nspintrct nspproteinrep nspbackbone nspdstl nspdataio nspgeometry nsploopclosure nsppdbstatistics nspdesignseq)
set(ALL_LIBS -Wl,--start-group ${SCUBA_LIBS} -Wl,--end-group ${PythonLibs_LIBRARIES} ${LAPACK_LIBRARIES})

set(USR_SRC_DIR ${MY_ROOT}/pySCUBA)
add_subdirectory (3rdparty/pybind11)
add_subdirectory (pySCUBA/dstl)
add_subdirectory (pySCUBA/dataio)
add_subdirectory (pySCUBA/geometry)
add_subdirectory (pySCUBA/loopclosure)
add_subdirectory (pySCUBA/backbone)
add_subdirectory (pySCUBA/pdbstatistics)
add_subdirectory (pySCUBA/proteinrep)
add_subdirectory (pySCUBA/designseq)
add_subdirectory (pySCUBA/iblock)
add_subdirectory (pySCUBA/sampling)
add_subdirectory (pySCUBA/pybind11)

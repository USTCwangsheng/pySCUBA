/*
 * @Descripttion: 
 * @version: 0.0
 * @Author: yuhaihan
 * @Date: 2022-04-06 09:48:16
 * @LastEditTime: 2022-04-09 20:56:40
 * fixdis_randomorientation
 * samplerelativeposition 未成功
 */
#include <pybind11/pybind11.h>

#include "geometry/relativeposition.h"

namespace py = pybind11;
namespace geo=NSPgeometry;

PYBIND11_MODULE(pb_geometry, m) {
    py::class_<geo::RelativePosition>(m, "RelativePosition")
        .def(py::init<>())
        .def_readwrite("location",&geo::RelativePosition::location)
        .def_readwrite("orientation",&geo::RelativePosition::orientation)
        .def("read", &geo::RelativePosition::read)
        .def("write", &geo::RelativePosition::write)
        .def("tovector", &geo::RelativePosition::tovector);

    // m.def("fixdis_randomorientation",&geo::fixdis_randomorientation<double &>);
    // m.def("samplerelativeposition",&geo::samplerelativeposition<double &>);

};
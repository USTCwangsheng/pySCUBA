/*
 * @Descripttion: 
 * @version: 0.0
 * @Author: yuhaihan
 * @Date: 2022-04-06 09:41:20
 * @LastEditTime: 2022-04-09 19:58:34
 */
/*
 * @Descripttion: 
 * @version: 0.0
 * @Author: yuhaihan
 * @Date: 2022-04-06 09:41:20
 * @LastEditTime: 2022-04-06 09:41:21
 */
#include <pybind11/pybind11.h>

#include "geometry/rectgrid.h"

namespace py = pybind11;
namespace geo=NSPgeometry;

PYBIND11_MODULE(pb_geometry, m) {
    py::class_<geo::RectGrid>(m, "RectGrid")
        .def(py::init<const geo::XYZ &,const std::array<int,3> &,double>())
        .def(py::init<const geo::XYZ &,const geo::XYZ &,double>())
        .def_readwrite("p000",&geo::RectGrid::p000)
        .def_readwrite("step",&geo::RectGrid::step)
        .def_readwrite("sizes",&geo::RectGrid::sizes)
        .def("center", &geo::RectGrid::center)
        .def("idx2xyz", &geo::RectGrid::idx2xyz)
        .def("mindis2", &geo::RectGrid::mindis2)
        .def("neighborgridpoints", &geo::RectGrid::neighborgridpoints)
        .def("neighborxyzs", &geo::RectGrid::neighborxyzs);
};
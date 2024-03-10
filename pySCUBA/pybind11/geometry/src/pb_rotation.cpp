/*
 * @Descripttion: 
 * @version: 0.0
 * @Author: yuhaihan
 * @Date: 2022-04-06 10:11:13
 * @LastEditTime: 2022-04-08 14:29:42
 * init函数没有必要重载，每一个都在初始化时自动调用
 * using lambda function to center
 */
#include <pybind11/pybind11.h>

#include "geometry/rotation.h"

namespace py = pybind11;
namespace geo=NSPgeometry;

PYBIND11_MODULE(pb_geometry, m) {
    py::class_<geo::Rotation>(m, "Rotation")
        .def(py::init<>())
        .def(py::init<const geo::XYZ &,const geo::XYZ &, const geo::XYZ & , const geo::XYZ &,double>())
        .def(py::init<const geo::QuaternionCrd &, const geo::XYZ &>())
        .def("apply",&geo::Rotation::apply)
        .def("center",
            [](geo::Rotation &a) {
                geo::XYZ center;
                center=a.center();
                return center;
            }
        )
        .def("applytoCopy",&geo::Rotation::applytoCopy);
};
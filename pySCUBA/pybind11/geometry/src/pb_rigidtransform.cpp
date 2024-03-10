/*
 * @Descripttion: 
 * @version: 0.0
 * @Author: yuhaihan
 * @Date: 2022-04-06 09:59:46
 * @LastEditTime: 2022-04-12 10:39:31
 */
#include <pybind11/pybind11.h>

#include "geometry/rigidtransform.h"
#include "geometry/rotation.h"

namespace py = pybind11;
namespace geo=NSPgeometry;

PYBIND11_MODULE(pb_geometry, m) {
    py::class_<geo::RigidTransform>(m, "RigidTransform")
        .def(py::init<const geo::QuaternionCrd &, const geo::XYZ &, const geo::XYZ &>())
        .def(py::init<const geo::Rotation &>())
        .def(py::init<const geo::Rotation &, const geo::XYZ &>())
        .def(py::init<const geo::XYZ & , const geo::XYZ &, const geo::XYZ &, const geo::XYZ &,double , double , double >())
        .def(py::init<const std::vector<geo::XYZ> &, const std::vector<double> &>())
        .def("axistoorigin", &geo::RigidTransform::axistoorigin)
        .def("getreverse", &geo::RigidTransform::getreverse)
        .def("apply", &geo::RigidTransform::apply)
        .def("translation",
            [](geo::RigidTransform &a) {
                geo::XYZ translation=a.translation();
                return translation;
            }
        )
        .def("rotation",
            [](geo::RigidTransform &a) {
                return a.rotation();
            }
        )
        .def("applytoCopy", &geo::RigidTransform::applytoCopy);
};
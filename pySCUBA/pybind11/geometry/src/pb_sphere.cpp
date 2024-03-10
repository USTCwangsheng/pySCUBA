/*
 * @Descripttion: 
 * @version: 0.0
 * @Author: yuhaihan
 * @Date: 2022-04-06 10:18:36
 * @LastEditTime: 2022-04-08 23:30:37
 * randompoint 模板 里的RNG模板类型使用double，该处理待保留意见
 */
#include <pybind11/pybind11.h>

#include "geometry/sphere.h"

namespace py = pybind11;
namespace geo=NSPgeometry;

PYBIND11_MODULE(pb_geometry, m) {
    py::class_<geo::Sphere>(m, "Sphere")
        .def(py::init<>())
        .def(py::init<const geo::XYZ &, double >())
        .def_readwrite("center",&geo::Sphere::center)
        .def_readwrite("radius2",&geo::Sphere::radius2)
        .def("randompoint",&geo::Sphere::randompoint<double &>)
        .def("insphere", &geo::Sphere::insphere);
};
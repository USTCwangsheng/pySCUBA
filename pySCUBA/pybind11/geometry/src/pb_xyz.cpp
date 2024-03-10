/*
 * @Descripttion: 
 * @version: 0.0
 * @Author: yuhaihan
 * @Date: 2022-04-02 16:07:58
 * @LastEditTime: 2022-04-12 11:23:39
 */
/*
 * @Descripttion: 
 * @version: 0.0
 * @Author: yuhaihan
 * @Date: 2022-03-26 21:21:55
 * @LastEditTime: 2022-03-30 17:02:45
 */

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include "geometry/xyz.h"

namespace py = pybind11;
namespace geo=NSPgeometry;

PYBIND11_MODULE(pb_geometry, m) {
    py::class_<geo::XYZ>(m, "XYZ")
        .def(py::init<double, double, double>(), py::arg("ix") = 0.0, py::arg("iy") = 0.0,py::arg("iz")=0.0)
        .def(py::init<const std::vector<double> &>())
        .def(py::init<double &,double>())
        .def_readwrite("x_",&geo::XYZ::x_)
        .def_readwrite("y_",&geo::XYZ::y_)
        .def_readwrite("z_",&geo::XYZ::z_)
        .def("squarednorm", &geo::XYZ::squarednorm)
        .def("length", &geo::XYZ::length)
        .def("squaredDistance", &geo::XYZ::squaredDistance)
        .def("distance", &geo::XYZ::distance)
        .def("__repr__", 
            [](geo::XYZ &a) {
                const std::string FMT {"%8.3f%8.3f%8.3f"};
	    	    char s[30];
		        sprintf(s, FMT.c_str(), a.x_, a.y_, a.z_);
                return "XYZ ("+std::string(s) + ")";
            }
        )
        .def(py::self + py::self)
        .def(py::self - py::self)
        // .def(py::self += py::self)
        .def(double() * py::self)
        .def(py::self * double())
        .def(py::self / double())
        // .def(py::self * float())
        .def("dot", &geo::dot)
        .def("cross", &geo::cross)
        .def(-py::self);
};
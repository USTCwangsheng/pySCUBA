/*
 * @Descripttion: 
 * @version: 0.0
 * @Author: yuhaihan
 * @Date: 2022-03-30 23:08:53
 * @LastEditTime: 2022-04-08 21:22:04
 */

#include <pybind11/pybind11.h>

#include "geometry/localframe.h"

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

namespace py = pybind11;
namespace geo=NSPgeometry;

PYBIND11_MODULE(pb_geometry, m) {
    py::class_<geo::LocalFrame>(m, "LocalFrame")
        .def(py::init<>())
        .def("local2globalcrd", &geo::LocalFrame::local2globalcrd)
        .def("global2localcrd",&geo::LocalFrame::global2localcrd)
        .def("axisline",&geo::LocalFrame::axisline)
        .def("make_localframe", overload_cast_<geo::XYZ,geo::XYZ,geo::XYZ>()(&geo::make_localframe))
        .def("make_localframe", overload_cast_<const geo::Line &, const geo::XYZ & >()(&geo::make_localframe))
        .def("make_localframe", overload_cast_<const geo::QuaternionCrd & ,const geo::XYZ &>()(&geo::make_localframe),py::arg("ori")=geo::XYZ(.0,.0,.0))
        .def_readwrite("origin",&geo::LocalFrame::origin_);


    py::class_<geo::A3LocalFrame,geo::LocalFrame>(m, "A3LocalFrame")
        .def(py::init<>())
        .def_readwrite("dlfda",&geo::A3LocalFrame::dlfda)
        .def("distributedvdx",&geo::A3LocalFrame::distributedvdx);

    py::class_<geo::DLocalFrameDx>(m, "DLocalFrameDx")
        .def(py::init<>())
        .def_readwrite("doridx",&geo::DLocalFrameDx::doridx)
        .def_readwrite("daxisdx",&geo::DLocalFrameDx::daxisdx);
};

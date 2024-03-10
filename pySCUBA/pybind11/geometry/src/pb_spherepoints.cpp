/*
 * @Descripttion: 
 * @version: 0.0
 * @Author: yuhaihan
 * @Date: 2022-04-06 10:23:36
 * @LastEditTime: 2022-04-09 21:37:37
 * 真是没想到，使用py::arg("")=,只检查数目，不检查名称，这真的可以成功初始化吗
 */
#include <pybind11/pybind11.h>

#include "geometry/spherepoints.h"


template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

namespace py = pybind11;
namespace geo=NSPgeometry;

PYBIND11_MODULE(pb_geometry, m) {
    m.def("genspherepoints",overload_cast_<int, std::vector<std::vector<double>> &,double >()(&geo::genspherepoints),py::arg("a")=0,py::arg("b")=0,py::arg("r")=1.0);
    m.def("genspherepoints",overload_cast_<int, std::vector<geo::XYZ> &,double >()(&geo::genspherepoints),py::arg("a")=0,py::arg("b")=0,py::arg("r")=1.0);
    m.def("genspherepoints",overload_cast_<const geo::LocalFrame &,int,std::vector<geo::XYZ> &,double >()(&geo::genspherepoints),py::arg("a")=0,py::arg("c")=0,py::arg("b")=0,py::arg("r")=1.0);
};
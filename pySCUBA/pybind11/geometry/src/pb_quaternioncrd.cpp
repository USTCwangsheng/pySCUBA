/*
 * @Descripttion: 
 * @version: 0.0
 * @Author: yuhaihan
 * @Date: 2022-04-02 17:06:56
 * @LastEditTime: 2022-04-08 22:50:42
 * QuaternionCrd(XYZ xyz,double angle,double deg=3.14159265358979323846/180.0)
 * 以上一行不指定前两个的初始值会报错
 * getmatrix非常顽固的不行
 * operator[] 使用python 替代对象 __getitem__ 配合lambda函数完成
 */
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include "geometry/quaternioncrd.h"
#include "geometry/localframe.h"

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

namespace py = pybind11;
namespace geo=NSPgeometry;

PYBIND11_MODULE(pb_geometry, m) {
    py::class_<geo::QuaternionCrd>(m, "QuaternionCrd")
        .def(py::init<>())
        .def(py::init<double,int>())
        .def(py::init<geo::LocalFrame &>())
        .def(py::init<const std::vector<double> &>())
        .def(py::init<geo::XYZ,double,double>(),py::arg("xyz")=geo::XYZ(0,0,0),py::arg("angle")=0,py::arg("deg")=0.017453292519943295)
        .def("scale_rotation",&geo::QuaternionCrd::scale_rotation)
        .def("angle",&geo::QuaternionCrd::angle)
        .def("maptoR3",&geo::QuaternionCrd::maptoR3)
        .def("angle_rad",&geo::QuaternionCrd::angle_rad)
        .def("axis",&geo::QuaternionCrd::axis)
        .def("vectcomp",&geo::QuaternionCrd::vectcomp)
        // .def("getmatrix",overload_cast_<double (&)[3][3]>()(&geo::QuaternionCrd::getmatrix, py::const_))
        .def("diff",overload_cast_<geo::QuaternionCrd>()(&geo::QuaternionCrd::diff))
        .def("diff",overload_cast_<std::vector<double>,std::vector<double>>()(&geo::QuaternionCrd::diff))
        .def("shift",&geo::QuaternionCrd::shift)
        .def("mean",&geo::QuaternionCrd::mean)
        .def("invert",&geo::QuaternionCrd::invert)
        .def("negative_equivalent",&geo::QuaternionCrd::negative_equivalent)
        .def("diff_rad",&geo::QuaternionCrd::diff_rad)
        .def(py::self * py::self)
        .def("quaternionaligntwovectors",&geo::quaternionaligntwovectors)
        .def("__getitem__",
            [](geo::QuaternionCrd &a,int i) {
                return a.Q[i];
            }
        )
        .def_readwrite("Q", &geo::QuaternionCrd::Q);
};
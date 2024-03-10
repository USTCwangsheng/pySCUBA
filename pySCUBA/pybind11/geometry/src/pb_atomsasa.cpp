/*
 * @Descripttion: 
 * @version: 0.0
 * @Author: yuhaihan
 * @Date: 2022-03-30 14:18:09
 * @LastEditTime: 2022-04-08 16:12:29
 * init 应该不需要定义接口
 * crd函数使用lambda作为替代
 */
#include <pybind11/pybind11.h>

#include "geometry/localframe.h"
#include "geometry/atomsasa.h"

namespace py = pybind11;
namespace geo=NSPgeometry;

PYBIND11_MODULE(pb_geometry, m) {
    py::class_<geo::SASAParameter>(m, "SASAParameter")
        .def(py::init<double,double,int>(), py::arg("solventradius") = 1.4, py::arg("soluteradius") = 1.4,py::arg("nsurfacepoints")=256)
        .def_readwrite("solventradius", &geo::SASAParameter::solventradius)
        .def_readwrite("soluteradius",&geo::SASAParameter::soluteradius)
        .def_readwrite("nsurfacepoints",&geo::SASAParameter::nsurfacepoints);
    py::class_<geo::AtomSASA>(m, "AtomSASA")
        .def(py::init<>())
        .def(py::init<const geo::SASAParameter &>())
        .def(py::init<const geo::LocalFrame &, const geo::SASAParameter &>())
        .def("crd",
            [](geo::AtomSASA &a) {
                geo::XYZ crd=a.crd();
                return crd;
            }
        )
        .def("par",
            [](geo::AtomSASA &a) {
                geo::SASAParameter par=a.par();
                return par;
            }
        )
        .def("update",&geo::AtomSASA::update)
        .def("updateAtomSASA",&geo::updateAtomSASA);
};

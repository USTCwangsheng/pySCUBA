/*
 * @Descripttion: 
 * @version: 0.0
 * @Author: yuhaihan
 * @Date: 2022-04-04 20:25:35
 * @LastEditTime: 2022-04-04 20:57:14
 */
#include <pybind11/pybind11.h>

#include "sampling/groupdistance.h"
#include "geometry/calculators.h"

namespace py = pybind11;
namespace spl=NSPsampling;

PYBIND11_MODULE(pb_sampling, m) {
    py::class_<spl::GroupContact>(m, "GroupContact")
        // .def("theta",static_cast<double (spl::GroupContact::*)(const std::vector<NSPgeometry::XYZ> &,std::vector<DvDxi> *)>(&spl::GroupContact::theta))
        .def("theta",static_cast<double const(spl::GroupContact::*)(const double,const double *)>(&spl::GroupContact::theta));
};
/*
 * @Descripttion: 
 * @version: 0.0
 * @Author: yuhaihan
 * @Date: 2022-04-05 23:13:24
 * @LastEditTime: 2022-04-06 09:40:20
 * setup,fitting 使用lambda函数绕过
 */
#include <pybind11/pybind11.h>

#include "geometry/quatfit.h"

namespace py = pybind11;
namespace geo=NSPgeometry;

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;


PYBIND11_MODULE(pb_geometry, m) {
    py::class_<geo::QuatFit>(m, "QuatFit")
        .def(py::init<>())
        // .def("setup",overload_cast_<const std::vector<double> &,const std::vector<double> &,std::vector<double>>()(&geo::QuatFit::setup))
        // .def("setup",overload_cast_<std::vector<double> &,std::vector<double> &,const std::vector<double>>()(&geo::QuatFit::setup))
        .def("setup",
            [](geo::QuatFit &a,std::vector<double> &ref_coord,std::vector<double> & coord,std::vector<double> weights=std::vector<double>()){
                return a.setup(ref_coord,coord,weights);
            })
        .def("evec",&geo::QuatFit::evec)
        .def("matrix",&geo::QuatFit::matrix)
        .def("getquaternioncrd",&geo::QuatFit::getquaternioncrd)
        // .def("fitting",overload_cast_<const std::vector<double> &,const std::vector<double> &,std::vector<double>>()(&geo::QuatFit::fitting))
        // .def("fitting",overload_cast_<const std::vector<geo::XYZ> &,const std::vector<geo::XYZ> &,std::vector<double>>()(&geo::QuatFit::fitting))
        .def("fitting",
            [](geo::QuatFit &a,std::vector<double> &ref_coord,std::vector<double> & coord,std::vector<double> weights=std::vector<double>()){
                return a.fitting(ref_coord,coord,weights);
            })
        .def("transform",&geo::QuatFit::transform)
        .def("RMSD",&geo::RMSD)
        .def("getRigidTransform",&geo::QuatFit::getRigidTransform);
};
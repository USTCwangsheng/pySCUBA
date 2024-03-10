/*
 * @Descripttion: 
 * @version: 0.0
 * @Author: yuhaihan
 * @Date: 2022-04-01 15:38:39
 * @LastEditTime: 2022-04-08 23:20:56
 * meansquareddist2 未解决 因为模板
 */
#include <pybind11/pybind11.h>

#include "geometry/calculators.h"
#include "geometry/xyz.h"


namespace py = pybind11;
namespace geo=NSPgeometry;

PYBIND11_MODULE(pb_geometry, m) {
    m.def("distance", static_cast<double (*)(const geo::XYZ &,const geo::XYZ &)>(&geo::distance));
    m.def("distance", static_cast<double (*)(const geo::XYZ &,const geo::XYZ &,std::vector<geo::XYZ> *)>(&geo::distance));
    m.def("distance2", &geo::distance2);
    m.def("angle", static_cast<double (*)(const geo::XYZ &,const geo::XYZ &,const geo::XYZ &)>(&geo::angle));    
    m.def("angle", static_cast<double (*)(const geo::XYZ &,const geo::XYZ &,const geo::XYZ &,std::vector<geo::XYZ> *)>(&geo::angle));  
    m.def("cos_angle", &geo::cos_angle);
    m.def("torsion", static_cast<double (*)(const geo::XYZ &,const geo::XYZ &,const geo::XYZ &,const geo::XYZ &,std::vector<geo::XYZ> *)>(&geo::torsion));    
    m.def("torsion", static_cast<double (*)(const geo::XYZ &,const geo::XYZ &,const geo::XYZ &,const geo::XYZ &,std::vector<geo::XYZ> *)>(&geo::torsion));    
    m.def("rmsd", &geo::rmsd);
    m.def("center",&geo::center);
    m.def("InternaltoXYZ", static_cast<geo::XYZ (*)(const geo::XYZ &,const geo::XYZ &,const geo::XYZ &,double,double,double)>(&geo::InternaltoXYZ));
    m.def("InternaltoXYZ", static_cast<geo::XYZ (*)(const geo::XYZ &,const geo::XYZ &,double,double)>(&geo::InternaltoXYZ));
    m.def("InternaltoXYZ", static_cast<geo::XYZ (*)(const geo::XYZ &,double)>(&geo::InternaltoXYZ));
    m.def("radiusgyr", static_cast<double (*)(const std::vector<NSPgeometry::XYZ> & )>(&geo::radiusgyr));
    m.def("radiusgyr", static_cast<double (*)(const std::vector<double> &)>(&geo::radiusgyr));
};
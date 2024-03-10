/*
 * @Descripttion: 
 * @version: 0.0
 * @Author: yuhaihan
 * @Date: 2022-04-06 10:32:58
 * @LastEditTime: 2022-04-06 10:53:25
 */
#include <pybind11/pybind11.h>

#include "geometry/structalign.h"

namespace py = pybind11;
namespace geo=NSPgeometry;

PYBIND11_MODULE(pb_geometry, m) {
    m.def("alignpointset", &geo::alignpointset);
};
/*
 * @Descripttion: 
 * @version: 0.0
 * @Author: yuhaihan
 * @Date: 2022-04-02 16:02:14
 * @LastEditTime: 2022-04-02 16:07:33
 */
#include <pybind11/pybind11.h>


#include "geometry/loopccd.h"

namespace py = pybind11;
namespace geo=NSPgeometry;

PYBIND11_MODULE(pb_geometry, m) {
    m.def("ccdangle", &geo::ccdangle);
    m.def("rotationtormsd", &geo::rotationtormsd);
};
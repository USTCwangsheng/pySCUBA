#include "backbone/backbonebuilder.h"
#include <pybind11/pybind11.h>

namespace py=pybind11;
using namespace NSPproteinrep;

PYBIND11_MODULE(pybackbonebuilder,m)
{
    py::class_<BackBoneBuilder>(m, "BackBoneBuilder")
        // .def(py::init<double &>())
        .def("buildforwardbackbone",&BackBoneBuilder::buildforwardbackbone);
}
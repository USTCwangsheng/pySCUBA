#include "backbone/backbonesite.h"
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace std;
using namespace NSPproteinrep;


PYBIND11_MODULE(pybackbonesite,m)
{
    py::class_<BackBoneSite>(m, "BackBoneSite", pybind11::buffer_protocol())
        // .def(py::init<string, string, char, char, int, int, bool, vector<double>>());
        .def_readwrite("resname", &BackBoneSite::resname)
        .def_readwrite("pdbid",&BackBoneSite::pdbid)
        .def_readwrite("chainid",&BackBoneSite::chainid)
        .def_readwrite("sscode",&BackBoneSite::sscode)
        .def_readwrite("resid",&BackBoneSite::resid)
        .def_readwrite("resseq",&BackBoneSite::resseq)
        .def_readwrite("isgap",&BackBoneSite::isgap)
        // .def_readwrite("data_",&BackBoneSite::data_)
        .def("sscodechar",&BackBoneSite::sscodechar);
}
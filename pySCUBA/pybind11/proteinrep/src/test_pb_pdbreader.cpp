/*
 * Author: cyx
 * 2022-03-31
 */

#include <pybind11/pybind11.h>
#include "proteinrep/pdbreader.h"
namespace prp = NSPproteinrep;
namespace py = pybind11;

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

PYBIND11_MODULE(pb_proteinrep, m)
{
	py::class_<prp::PdbReader>(m, "PdbReader")
//	    .def("readpdb", overload_cast_<const std::string &>()(&prp::PdbReader::readpdb))
//	    .def("readpdb", overload_cast_<std::vector<std::string> &>()(&prp::PdbReader::readpdb))
//	    .def("addRecord", &prp::PdbReader::addRecord);
}

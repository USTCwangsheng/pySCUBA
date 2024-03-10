/*
 * @Author: cyx
 * @2022-03-31
 */

#include <pybind11/pybind11.h>
#include "proteinrep/aaconformer.h"
namespace py = pybind11;
namespace prp = NSPproteinrep;

/*
template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

prp::AAConformer ma(const std::vector<prp::PdbRecord> & residuerecords)
{
	return make_aaconformer(residuerecords);
}
*/

PYBIND11_MODULE(pb_proteinrep, m)
{
	py::class_<prp::AAConformer>(m, "AAConformer")
//	    .def("make_aaconformer", &ma)
//	    .def_readwrite_static("sidechainatoms", &prp::AAConformer::sidechainatoms)
//	    .def_readwrite_static("backboneatoms", &prp::AAConformer::backboneatoms)
//	    .def_readwrite_static("mainchainatoms", &prp::AAConformer::mainchainatoms)
	    .def_readwrite("residuename", &prp::AAConformer::residuename)
	    .def_readwrite("chainid_or", &prp::AAConformer::chainid_or)
	    .def_readwrite("residueid_or", &prp::AAConformer::residueid_or)
	    .def_readwrite("atomlist", &prp::AAConformer::atomlist)
	    .def_readwrite("globalcrd", &prp::AAConformer::globalcrd)
	    .def_readwrite("localcrd", &prp::AAConformer::localcrd)
	    .def_readwrite("localframe", &prp::AAConformer::localframe);

	py::class_<prp::AAConformersInModel>(m, "AAConformersInModel")
//	    .def("readpdbfile", &prp::AAConformersInModel::readpdbfile)
//	    .def("getconformers", &prp::AAConformersInModel::getconformers)
	    .def_readwrite("conformers", &prp::AAConformersInModel::conformers);

}

#include "dstl/graph.h"
#include "dstl/globalinstances.h"
#include "dstl/tree.h"
#include "dstl/stlutil.h"
#include "dstl/domaintree.h"
#include "dstl/ga.h"
#include "dstl/minimizer.h"
#include "dstl/nbstree.h"
#include "dstl/topn.h"
#include "dstl/randomengine.h"
#include "dstl/vec2d.h"
#include "dstl/permutation.h"
#include "dstl/sortindex.h"
#include "dstl/nestedcontainers.h"
#include "dstl/perwalker.h"
#include "dstl/symmatrix1d.h"
#include "dstl/paretofront.h"
#include "dstl/domainleaf.h"
#include "dstl/pca.h"
#include "dstl/nrtopn.h"
#include "dstl/alignset.h"
#include "dstl/vectortree.h"
#include "dstl/mapkeyint.h"
#include "dstl/nnearest.h"
#include "dstl/nnregressionmodel.h"

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;
namespace py = pybind11;
namespace geo=NSPgeometry;
namespace dstl=NSPdstl;

void reseed(int seed){
    NSPdstl::RandomEngine<>::getinstance().reseed(seed);
}

PYBIND11_MODULE(pb_dstl, m) {
    py::class_<dstl::RandomEngine<int>>(m,"RandomEngine")
        .def(py::init<>())
        .def("reseed",&reseed);
}




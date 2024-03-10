/*
 * @Descripttion: 
 * @version: 0.0
 * @Author: yuhaihan
 * @Date: 2022-04-01 21:27:14
 * @LastEditTime: 2022-04-12 14:34:26
 * 少量被注释掉的没有完成
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>

#include "geometry/crdtree.h"
#include "dstl/tree.h"
#include "geometry/calculators.h"

namespace py = pybind11;
namespace geo=NSPgeometry;


PYBIND11_MODULE(pb_geometry, m) {
    struct AIJK;
    typedef std::vector<std::string> ATOMKEY;
    typedef NSPdstl::TreeNode<ATOMKEY,AIJK> TreeNode;
    py::class_<geo::TopoTree<ATOMKEY>>topotree(m,"TopoTree");
    topotree.def(py::init<const ATOMKEY & >())
        // .def("attachAtom",&geo::TopoTree<ATOMKEY &>::attachAtom<const ATOMKEY &,const ATOMKEY &>)
        // .def("attachBranch",&geo::TopoTree<ATOMKEY &>::attachBranch)
        .def("tree",
            [](geo::TopoTree<ATOMKEY> &a) {
                return a.tree();
            }
        )
        .def("subtreemap",
            [](geo::TopoTree<ATOMKEY> &a) {
                return a.subtreemap();
            }
        )
        .def("treemapped",
            [](geo::TopoTree<ATOMKEY> &a) {
                return a.treemapped();
            }
        )
        .def("branch",
            [](geo::TopoTree<ATOMKEY> &a,ATOMKEY key) {
                return a.branch(key);
            }
        );

    py::class_<geo::TopoTree<ATOMKEY>::AIJK>(topotree, "AIJK")
        .def(py::init<const ATOMKEY &,const ATOMKEY &, const ATOMKEY &>())
        .def(py::init<const ATOMKEY &>())
        .def_readwrite("ia",&geo::TopoTree<ATOMKEY>::AIJK::ia)
        .def_readwrite("ja",&geo::TopoTree<ATOMKEY>::AIJK::ja)
        .def_readwrite("ka",&geo::TopoTree<ATOMKEY>::AIJK::ka)
        .def("appropriate",&geo::TopoTree<ATOMKEY>::AIJK::appropriate)
        .def("setParent",&geo::TopoTree<ATOMKEY>::AIJK::setParent);

    py::class_<geo::IntCrd>(m,"IntCrd")
        .def(py::init<double,double,double>(),py::arg("d")=-1.,py::arg("a")=10e3,py::arg("t")=10e3)
        .def(py::init<const geo::TopoTree<ATOMKEY> * ,const std::map<ATOMKEY,geo::XYZ> & ,const ATOMKEY &>())
        .def_readwrite("distance",&geo::IntCrd::distance)
        .def_readwrite("angle",&geo::IntCrd::angle)
        .def_readwrite("torsion",&geo::IntCrd::torsion)
        .def("validAngle",&geo::IntCrd::validAngle)
        .def("validTorsion",&geo::IntCrd::validTorsion)
        .def("validDistance",&geo::IntCrd::validDistance)
        .def("valid",&geo::IntCrd::valid);
    
    py::class_<geo::CrdTree<ATOMKEY>>(m,"CrdTree")
        .def(py::init<geo::TopoTree<ATOMKEY> *>())
        // .def("intcrdmap",&geo::CrdTree<ATOMKEY>::intcrdmap)
        .def("resetCrdMap",&geo::CrdTree<ATOMKEY>::resetCrdMap)
        .def("crdValid",&geo::CrdTree<ATOMKEY>::crdValid)
        .def("resetIntCrdMap",&geo::CrdTree<ATOMKEY>::resetIntCrdMap)
        .def("resetMaps",&geo::CrdTree<ATOMKEY>::resetMaps)
        .def("calcXYZ",&geo::CrdTree<ATOMKEY>::calcXYZ)
        .def("calcInternal",&geo::CrdTree<ATOMKEY>::calcInternal)
        .def("setDistance",&geo::CrdTree<ATOMKEY>::setDistance)
        .def("setAngle",&geo::CrdTree<ATOMKEY>::setAngle)
        .def("setTorsion",&geo::CrdTree<ATOMKEY>::setTorsion)
        .def("rotateBond",&geo::CrdTree<ATOMKEY>::rotateBond)
        // .def("setInternal",&geo::CrdTree<ATOMKEY>::setInternal)
        .def("updateInternal",&geo::CrdTree<ATOMKEY>::updateInternal)
        .def("addAtom",&geo::CrdTree<ATOMKEY>::addAtom)
        .def("independentXYZ",&geo::CrdTree<ATOMKEY>::independentXYZ)
        .def("copyCrdMap",&geo::CrdTree<ATOMKEY>::copyCrdMap)
        .def("startupdate",&geo::CrdTree<ATOMKEY>::startupdate)
        .def("intcrdmap",
            [](geo::CrdTree<ATOMKEY> &a) {
                return a.intcrdmap();
            }
        )
        .def("crdmap",
            [](geo::CrdTree<ATOMKEY> &a) {
                return a.crdmap();
            }
        )
        .def("topotree",
            [](geo::CrdTree<ATOMKEY> &a) {
                return a.topotree();
            }
        );

    typedef ATOMKEY KeyType;
    typedef typename geo::TopoTree<KeyType>::TreeNode Node;
    py::class_<geo::CalcNodeXYZ<ATOMKEY>>(m,"CalcNodeXYZ")
        .def(py::init<geo::CrdTree<KeyType> *,bool>())
        .def("enter",&geo::CalcNodeXYZ<ATOMKEY>::enter)
        .def("leave",&geo::CalcNodeXYZ<ATOMKEY>::leave);
};
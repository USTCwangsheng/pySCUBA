/*
 * @Descripttion: 
 * @version: 0.0
 * @Author: yuhaihan
 * @Date: 2022-03-26 21:21:55
 * @LastEditTime: 2022-04-12 14:24:01
 * radiusgyr undifuned symbol
 */

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>

#include "geometry/xyz.h"
#include "geometry/atomsasa.h"
#include "geometry/loopccd.h"
#include "geometry/localframe.h"
#include "geometry/calculators.h"
#include "geometry/quaternioncrd.h"
#include "geometry/structalign.h"
#include "geometry/sphere.h"
#include "geometry/spherepoints.h"
#include "geometry/quatfit.h"
#include "geometry/rectgrid.h"
#include "geometry/relativeposition.h"
#include "geometry/rigidtransform.h"
#include "geometry/rotation.h"
#include "geometry/crdtree.h"

#include "dstl/tree.h"

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;
namespace py = pybind11;
namespace geo=NSPgeometry;

double rmsd(std::vector<geo::XYZ> crds1, std::vector<geo::XYZ> crds2)
{
	geo::QuatFit qf;
	double r2 = qf.setup(crds1, crds2);
	return sqrt(r2);
}


PYBIND11_MODULE(pb_geometry, m) {
    py::class_<geo::XYZ>(m, "XYZ")
        .def(py::init<>())
	.def(py::init<double, double, double>(), py::arg("ix") = 0.0, py::arg("iy") = 0.0,py::arg("iz")=0.0)
        .def(py::init<const std::vector<double> &>())
        .def(py::init<double &,double>())
        .def_readwrite("x",&geo::XYZ::x_)
        .def_readwrite("y",&geo::XYZ::y_)
        .def_readwrite("z",&geo::XYZ::z_)
        .def("squarednorm", &geo::XYZ::squarednorm)
        .def("length", &geo::XYZ::length)
        .def("squaredDistance", &geo::XYZ::squaredDistance)
        .def("distance", &geo::XYZ::distance)
        .def("__repr__", 
            [](geo::XYZ &a) {
                const std::string FMT {"%8.3f%8.3f%8.3f"};
	    	    char s[30];
		        sprintf(s, FMT.c_str(), a.x_, a.y_, a.z_);
                return "XYZ ("+std::string(s) + ")";
            }
        )
        .def(py::self + py::self)
        .def(py::self - py::self)
        // .def(py::self += py::self)
        .def(double() * py::self)
        .def(py::self * double())
        .def(py::self / double())
        // .def(py::self * float())
        .def("dot", &geo::dot)
        .def("cross", &geo::cross)
        .def(-py::self);

    // loopccd.h
    m.def("ccdangle", &geo::ccdangle);
    m.def("rotationtormsd", &geo::rotationtormsd);

    // atomsasa.h
    py::class_<geo::SASAParameter>(m, "SASAParameter")
        .def(py::init<double,double,int>(), py::arg("solventradius") = 1.4, py::arg("soluteradius") = 1.4,py::arg("nsurfacepoints")=256)
        .def_readwrite("solventradius", &geo::SASAParameter::solventradius)
        .def_readwrite("soluteradius",&geo::SASAParameter::soluteradius)
        .def_readwrite("nsurfacepoints",&geo::SASAParameter::nsurfacepoints);
    py::class_<geo::AtomSASA>(m, "AtomSASA")
        .def(py::init<>())
        .def(py::init<const geo::SASAParameter &>())
        .def(py::init<const geo::LocalFrame &, const geo::SASAParameter &>())
        .def("crd",
            [](geo::AtomSASA &a) {
                geo::XYZ crd=a.crd();
                return crd;
            }
        )
        .def("par",
            [](geo::AtomSASA &a) {
                geo::SASAParameter par=a.par();
                return par;
            }
        )
        .def("update",&geo::AtomSASA::update)
        .def("updateAtomSASA",&geo::updateAtomSASA);

    //localframe.h
    py::class_<geo::LocalFrame>(m, "LocalFrame")
        .def(py::init<>())
        .def("local2globalcrd", &geo::LocalFrame::local2globalcrd)
        .def("global2localcrd",&geo::LocalFrame::global2localcrd)
        .def("axisline",&geo::LocalFrame::axisline)
        .def("make_localframe", overload_cast_<geo::XYZ,geo::XYZ,geo::XYZ>()(&geo::make_localframe))
        .def("make_localframe", overload_cast_<const geo::Line &, const geo::XYZ & >()(&geo::make_localframe))
        .def("make_localframe", overload_cast_<const geo::QuaternionCrd & ,const geo::XYZ &>()(&geo::make_localframe),py::arg("ori")=geo::XYZ(.0,.0,.0))
        .def_readwrite("origin",&geo::LocalFrame::origin_);


    py::class_<geo::A3LocalFrame,geo::LocalFrame>(m, "A3LocalFrame")
        .def(py::init<>())
        .def_readwrite("dlfda",&geo::A3LocalFrame::dlfda)
        .def("distributedvdx",&geo::A3LocalFrame::distributedvdx);

    py::class_<geo::DLocalFrameDx>(m, "DLocalFrameDx")
        .def(py::init<>())
        .def_readwrite("doridx",&geo::DLocalFrameDx::doridx)
        .def_readwrite("daxisdx",&geo::DLocalFrameDx::daxisdx);

    //calculator.h
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

    //"geometry/quaternioncrd.h"
    py::class_<geo::QuaternionCrd>(m, "QuaternionCrd")
        .def(py::init<>())
        .def(py::init<double,int>())
        .def(py::init<geo::LocalFrame &>())
        .def(py::init<const std::vector<double> &>())
        .def(py::init<geo::XYZ,double,double>(),py::arg("xyz")=geo::XYZ(0,0,0),py::arg("angle")=0,py::arg("deg")=0.017453292519943295)
        .def("scale_rotation",&geo::QuaternionCrd::scale_rotation)
        .def("angle",&geo::QuaternionCrd::angle)
        .def("maptoR3",&geo::QuaternionCrd::maptoR3)
        .def("angle_rad",&geo::QuaternionCrd::angle_rad)
        .def("axis",&geo::QuaternionCrd::axis)
        .def("vectcomp",&geo::QuaternionCrd::vectcomp)
        // .def("getmatrix",overload_cast_<double (&)[3][3]>()(&geo::QuaternionCrd::getmatrix, py::const_))
        .def("diff",overload_cast_<geo::QuaternionCrd>()(&geo::QuaternionCrd::diff))
        .def("diff",overload_cast_<std::vector<double>,std::vector<double>>()(&geo::QuaternionCrd::diff))
        .def("shift",&geo::QuaternionCrd::shift)
        .def("mean",&geo::QuaternionCrd::mean)
        .def("invert",&geo::QuaternionCrd::invert)
        .def("negative_equivalent",&geo::QuaternionCrd::negative_equivalent)
        .def("diff_rad",&geo::QuaternionCrd::diff_rad)
        .def(py::self * py::self)
        .def("quaternionaligntwovectors",&geo::quaternionaligntwovectors)
        .def("__getitem__",
            [](geo::QuaternionCrd &a,int i) {
                return a.Q[i];
            }
        )
        .def_readwrite("Q", &geo::QuaternionCrd::Q);
    // struture align
    m.def("alignpointset", &geo::alignpointset);

    // sphere
    py::class_<geo::Sphere>(m, "Sphere")
        .def(py::init<>())
        .def(py::init<const geo::XYZ &, double >())
        .def_readwrite("center",&geo::Sphere::center)
        .def_readwrite("radius2",&geo::Sphere::radius2)
        .def("randompoint",&geo::Sphere::randompoint<double &>)
        .def("insphere", &geo::Sphere::insphere);

    // quatfit
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
    // rectgrid
    py::class_<geo::RectGrid>(m, "RectGrid")
        .def(py::init<const geo::XYZ &,const std::array<int,3> &,double>())
        .def(py::init<const geo::XYZ &,const geo::XYZ &,double>())
        .def_readwrite("p000",&geo::RectGrid::p000)
        .def_readwrite("step",&geo::RectGrid::step)
        .def_readwrite("sizes",&geo::RectGrid::sizes)
        .def("center", &geo::RectGrid::center)
        .def("idx2xyz", &geo::RectGrid::idx2xyz)
        .def("mindis2", &geo::RectGrid::mindis2)
        .def("neighborgridpoints", &geo::RectGrid::neighborgridpoints)
        .def("neighborxyzs", &geo::RectGrid::neighborxyzs);
    
    //relativeposition
    py::class_<geo::RelativePosition>(m, "RelativePosition")
        .def(py::init<>())
        .def_readwrite("location",&geo::RelativePosition::location)
        .def_readwrite("orientation",&geo::RelativePosition::orientation)
        .def("read", &geo::RelativePosition::read)
        .def("write", &geo::RelativePosition::write)
        .def("tovector", &geo::RelativePosition::tovector);

    // rigidtransform
    py::class_<geo::RigidTransform>(m, "RigidTransform")
        .def(py::init<const geo::QuaternionCrd &, const geo::XYZ &, const geo::XYZ &>())
        .def(py::init<const geo::Rotation &>())
        .def(py::init<const geo::Rotation &, const geo::XYZ &>())
        .def(py::init<const geo::XYZ & , const geo::XYZ &, const geo::XYZ &, const geo::XYZ &,double , double , double >())
        .def(py::init<const std::vector<geo::XYZ> &, const std::vector<double> &>())
        .def("axistoorigin", &geo::RigidTransform::axistoorigin)
        .def("getreverse", &geo::RigidTransform::getreverse)
        .def("apply", &geo::RigidTransform::apply)
        .def("translation",
            [](geo::RigidTransform &a) {
                geo::XYZ translation=a.translation();
                return translation;
            }
        )
        .def("rotation",
            [](geo::RigidTransform &a) {
                return a.rotation();
            }
        )
        .def("applytoCopy", &geo::RigidTransform::applytoCopy);

    //rotation
    py::class_<geo::Rotation>(m, "Rotation")
        .def(py::init<>())
        .def(py::init<const geo::XYZ &,const geo::XYZ &, const geo::XYZ & , const geo::XYZ &,double>())
        .def(py::init<const geo::QuaternionCrd &, const geo::XYZ &>())
        .def("apply",&geo::Rotation::apply)
        .def("center",
            [](geo::Rotation &a) {
                geo::XYZ center;
                center=a.center();
                return center;
            }
        )
        .def("applytoCopy",&geo::Rotation::applytoCopy);
    // spherepoints
    m.def("genspherepoints",overload_cast_<int, std::vector<std::vector<double>> &,double >()(&geo::genspherepoints),py::arg("a")=0,py::arg("b")=0,py::arg("r")=1.0);
    m.def("genspherepoints",overload_cast_<int, std::vector<geo::XYZ> &,double >()(&geo::genspherepoints),py::arg("a")=0,py::arg("b")=0,py::arg("r")=1.0);
    m.def("genspherepoints",overload_cast_<const geo::LocalFrame &,int,std::vector<geo::XYZ> &,double >()(&geo::genspherepoints),py::arg("a")=0,py::arg("c")=0,py::arg("b")=0,py::arg("r")=1.0);

    //crdtree.h
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
    typedef typename geo::TopoTree<KeyType> ::TreeNode Node;
    py::class_<geo::CalcNodeXYZ<ATOMKEY>>(m,"CalcNodeXYZ")
        .def(py::init<geo::CrdTree<KeyType> *,bool>())
        .def("enter",&geo::CalcNodeXYZ<ATOMKEY>::enter)
        .def("leave",&geo::CalcNodeXYZ<ATOMKEY>::leave);
};

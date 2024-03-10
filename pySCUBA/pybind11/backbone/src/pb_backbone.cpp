#include "backbone/backbonesite.h"
#include "backbone/backbonebuilder.h"
#include "backbone/backbonealignment.h"
// 导入pybind相关包
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
// #include <pybind11/complex.h>
// #include <pybind11/functional.h>
// #include <pybind11/chrono.h>

namespace py = pybind11;
using namespace std;
using namespace NSPproteinrep;

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

const double* get_data_(const BackBoneSite &bbs){   
    return bbs.data_;
}

void set_data_(BackBoneSite &bbs,const std::vector<double> data){
    assert(data.size()==BackBoneSite::DataField::DATADIM);
    for(int n=0;n<BackBoneSite::DataField::DATADIM;n++){
        bbs.data_[n]=data[n];
    }
}

typedef std::vector<std::pair<int,int>> AlignedPositions;

PYBIND11_MODULE(pb_backbone,m)
{
    // auto dtype = py::dtype(py::format_descriptor<double>::format());
    py::class_<BackBoneSite> backbonesite(m, "BackBoneSite");
    backbonesite.def_readwrite("resname", &BackBoneSite::resname)
        .def_readwrite("pdbid",&BackBoneSite::pdbid)
        .def_readwrite("chainid",&BackBoneSite::chainid)
        .def_readwrite("sscode",&BackBoneSite::sscode)
        .def_readwrite("resid",&BackBoneSite::resid)
        .def_readwrite("resseq",&BackBoneSite::resseq)
        .def_readwrite("isgap",&BackBoneSite::isgap)
        // .def_readwrite("data_",&BackBoneSite::data_)
        .def_property("data_",get_data_,set_data_)
        // .attr("data_") = py::array(dtype, {BackBoneSite::DataField::DATADIM}, {sizeof(double)}, array, nullptr)
        // .def("get_data", []() {
        //     return py::memoryview::from_buffer(
        //         &BackBoneSite::data_,                                    // buffer pointer
        //         {BackBoneSite::DataField::DATADIM},      // shape (rows, cols)
        //         { sizeof(double) * 4, sizeof(double) }   // strides in bytes
        //     );
        // })
        .def(py::init<>())

        //class variables
        .def("data", overload_cast_<>()(&BackBoneSite::data))
        .def("data", overload_cast_<>()(&BackBoneSite::data, py::const_))
        .def("phi", overload_cast_<>()(&BackBoneSite::phi, py::const_))
        .def("phi", overload_cast_<const BackBoneSite &>()(&BackBoneSite::phi))
        .def("psi", overload_cast_<>()(&BackBoneSite::psi, py::const_))
        .def("psi", overload_cast_<const BackBoneSite &>()(&BackBoneSite::psi))
        .def("omiga", overload_cast_<>()(&BackBoneSite::omiga, py::const_))
        .def("omiga", overload_cast_<const BackBoneSite &>()(&BackBoneSite::omiga))      
        
        //functions
        .def("nextpepcis",&BackBoneSite::nextpepcis)
        .def("settorsions",&BackBoneSite::settorsions)
        .def("sscodechar",&BackBoneSite::sscodechar)
        .def("__repr__", &BackBoneSite::toString)
        .def("localframe",&BackBoneSite::localframe)
        // .def("read",&BackBoneSite::read<const string::iterator&>)
        // template <typename... Ts>
        // template<typename IT>
        // .def("read", &BackBoneSite::read, py::call_guard<IT>());
        .def("getcrd", overload_cast_<int>()(&BackBoneSite::getcrd, py::const_))
        .def("getcrd", overload_cast_<std::vector<double> &>()(&BackBoneSite::getcrd, py::const_))
        .def("getcrd", overload_cast_<std::vector<NSPgeometry::XYZ> &>()(&BackBoneSite::getcrd, py::const_))
        .def("ncrd",&BackBoneSite::ncrd)
        .def("cacrd",&BackBoneSite::cacrd)
        .def("ccrd",&BackBoneSite::ccrd)
        .def("ocrd",&BackBoneSite::ocrd)
        
        // .def("changecrd", overload_cast_<const std::vector<double> &>()(&BackBoneSite::changecrd))
        // .def("changecrd", overload_cast_<const std::vector<NSPgeometry::XYZ> &>()(&BackBoneSite::changecrd))
        // .def("changecrd", overload_cast_<vector<NSPgeometry::XYZ>::iterator>()(&BackBoneSite::changecrd))
        .def("changecrd",[](BackBoneSite & bbs, const std::vector<double> &crd){
            bbs.changecrd(crd);
        }) //?
        
        .def("==", [](const BackBoneSite &) {
            return false;
            }, py::is_operator())

        .def("hcrd",&BackBoneSite::hcrd)
        .def("cd_procrd",&BackBoneSite::cd_procrd)
        .def("cbcrd",&BackBoneSite::cbcrd)
        .def("newocrdfrompsi",&BackBoneSite::newocrdfrompsi)
        .def("newncrdfrompsi",&BackBoneSite::newncrdfrompsi)
        .def("translate",&BackBoneSite::translate)
        .def("rotate",&BackBoneSite::rotate)
        .def("genPdbRecords",&BackBoneSite::genPdbRecords)
    ;

    py::enum_<BackBoneSite::DataField>(backbonesite, "Field")
        .value("PHI", BackBoneSite::DataField::PHI)
        .value("PSI", BackBoneSite::DataField::PSI)
        .value("OMIGA", BackBoneSite::DataField::OMIGA)
        .value("SASA", BackBoneSite::DataField::SASA)
        .value("NCRD", BackBoneSite::DataField::NCRD)
        .value("CACRD", BackBoneSite::DataField::CACRD)
        .value("CCRD", BackBoneSite::DataField::CCRD)
        .value("OCRD", BackBoneSite::DataField::OCRD)
        .value("DATADIM", BackBoneSite::DataField::DATADIM)
    ;

    m.def("readbackbonesites", overload_cast_<const std::string & ,std::vector<BackBoneSite> &>()(&readbackbonesites));
    // m.def("readbackbonesites", overload_cast_<std::istream &, int, std::vector<BackBoneSite> &>()(&readbackbonesites));
    m.def("writeSitesToPDB",[](std::string filename, std::vector<BackBoneSite> sites){
        std::ofstream os(filename);
        writeSitesToPDB(os, sites);
    });
    // m.def("readbackbonefrompdb", &readbackbonefrompdb);
    m.def("readbackbonefrompdb",[](std::string filename){
        std::vector<std::vector<BackBoneSite>> chains;
        readbackbonefrompdb(filename, chains);
        return chains;
    });
    m.def("generaterandombackbone", &generaterandombackbone);
    m.def("genbackbonesite", [](BackBoneSite psite, bool cispep,
		double phi, double psi, BackBoneSite &newsite){
        std::shared_ptr<BackBoneSite> newpsite=std::shared_ptr<BackBoneSite>(new BackBoneSite(psite));
        std::shared_ptr<BackBoneSite> newnewsite=std::shared_ptr<BackBoneSite>(new BackBoneSite(newsite));
        genbackbonesite(newpsite.get(), cispep, phi, psi, newnewsite.get());
        newsite=*(newnewsite.get());
    });
    m.def("genbackbonesite", [](bool cispep,
		double phi, double psi, BackBoneSite &newsite){
        std::shared_ptr<BackBoneSite> newnewsite=std::shared_ptr<BackBoneSite>(new BackBoneSite(newsite));
        genbackbonesite(nullptr, cispep, phi, psi, newnewsite.get());
        newsite=*(newnewsite.get());
    });
    m.def("genprevbackbonesite", &genprevbackbonesite);

    py::class_<BackBoneBuilder>(m, "BackBoneBuilder")
        .def(py::init<>())
        .def("buildstrandat",&BackBoneBuilder::buildstrandat)
        .def("buildhelixat",&BackBoneBuilder::buildhelixat)
        // .def("buildforwardbackbone",&BackBoneBuilder::buildforwardbackbone)
        .def("buildbackwardbackbone",&BackBoneBuilder::buildbackwardbackbone)
        // .def("buildlinkers",&BackBoneBuilder::buildlinkers)
        .def("buildlinkers",[](int length,
		const BackBoneSite &nflankingsite,
		const BackBoneSite &cflankingsite,
		const std::vector<std::pair<int,int>> & helixregions,
		const std::vector<std::pair<int,int>> & strandregions,
		const std::vector<int> & cissites){
            std::set<int> newcissites;
            for(auto n:cissites) newcissites.insert(n);
            BackBoneBuilder::buildlinkers(length, nflankingsite, cflankingsite, helixregions, strandregions, newcissites);
        })
        // .def("movechainto",&BackBoneBuilder::movechainto)
    ;

    m.def("buildforwardbackbone",[](int length,
		const BackBoneSite & nflankingsite, const std::vector<std::pair<int,int>> & helixregions,
		const std::vector<std::pair<int,int>> &strandregions, const std::vector<int> & cissites){
            std::set<int> newcissites;
            for(int n=0;n<cissites.size();n++) newcissites.insert(cissites[n]);
            return BackBoneBuilder::buildforwardbackbone(length, nflankingsite, helixregions,strandregions,newcissites);
    });

    // m.def("movechainto",&BackBoneBuilder::movechainto);
    m.def("movechainto",[](NSPgeometry::XYZ r0,NSPgeometry::XYZ direction,
		bool forward, std::vector<BackBoneSite> &chain){
            // std::vector<BackBoneSite> newchain=chain;
            BackBoneBuilder::movechainto(r0, direction, forward, chain);
            // for(int n=0;n<newchain.size();n++) std::cout<<newchain[n].cacrd().toString()<<std::endl;
            return chain;
    });
    

    py::class_<BackBoneAlignment> backbonealignment(m, "BackBoneAlignment");
        backbonealignment.def(py::init<int>())
        .def("init",&BackBoneAlignment::init)
        .def("align",[](const std::vector<BackBoneSite> &confa, const std::vector<BackBoneSite> &confb, int seedmode={BackBoneAlignment::AUTO_SEEDS}){
            BackBoneAlignment ba(seedmode);
            ba.init(confa,confb);
            std::shared_ptr<BackBoneAlignment::Results> optres=std::shared_ptr<BackBoneAlignment::Results> (new BackBoneAlignment::Results);
            ba.optalign(optres.get());
            return optres.get();
        })
        // .def("optalign",&BackBoneAlignment::optalign)
        // .def("selectseeds",&BackBoneAlignment::selectseeds)
        // .def("simple_selectseeds",&BackBoneAlignment::simple_selectseeds)
        // .def("SS_selectseeds",&BackBoneAlignment::SS_selectseeds)
        // .def("tryalign",&BackBoneAlignment::tryalign)
    ;

    // py::enum_<BackBoneAlignment>(backbonealignment, "")
    //     .value("AUTO_SEEDS", BackBoneAlignment::AUTO_SEEDS)
    //     .value("SIMPLE_SEEDS", BackBoneAlignment::SIMPLE_SEEDS)
    //     .value("SS_SEEDS", BackBoneAlignment::SS_SEEDS)        
    // ;

    py::class_<BackBoneAlignment::Results>(m, "Results")
        .def(py::init<>())
        .def_readonly("rigidtransform",&BackBoneAlignment::Results::rigidtransform)
        .def_readonly("alignedpositions",&BackBoneAlignment::Results::alignedpositions)
        // .def("length",&BackBoneAlignment::align)
        // .def("better",&BackBoneAlignment::optalign)
        
    ;
}
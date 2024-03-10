/*
 * Author: cyx
 * 2022-03-31
 */

#include <pybind11/pybind11.h>
#include "proteinrep/pdbrecord.h"
#include "proteinrep/intatomkey.h"
#include "dataio/inputlines.h"
#include "dataio/datapaths.h"
#include "geometry/xyz.h"
#include "proteinrep/aaconformer.h"
#include "proteinrep/intatomkey.h"
#include "proteinrep/residuestate.h"
#include "proteinrep/aminoacidseq.h"
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>

namespace py = pybind11;
using namespace std;
using namespace NSPproteinrep;
using namespace NSPgeometry;

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

//
PdbRecord mkrecord(unsigned int key, NSPgeometry::XYZ xyz,int atomid=0, int residueidshift=0,
		int elementsymbolsize=1){
	return make_pdbrecord<AtomKeyTypeA,NSPgeometry::XYZ>(key,xyz,atomid,residueidshift,elementsymbolsize);
}

//
std::string get_pdr_elementname(const PdbRecord &pdr){
	std::string res("  ");
	res[0]=pdr.elementname[0];
	res[1]=pdr.elementname[1];
	return res;
}
void set_pdr_elementname(PdbRecord &pdr,const std::string &name){
	        pdr.elementname[0]=name[0];
		        pdr.elementname[1]=name[1];
}

//
AAConformer pb_make_aaconformer(const std::vector<PdbRecord> & residuerecords, std::vector<AAConformer> *altconformers, bool replaceMSE)
{
	return make_aaconformer(residuerecords, altconformers, replaceMSE);
}
AAConformer pb_make_aaconformer(NSPdesignseq::Rotamer *rt, const BackBoneSite *site)
{
	return make_aaconformer(rt, site);
}

vector<XYZ> getcrdbyname(AAConformersInModel aim, string name)
{
	vector<XYZ> crds;
	for (auto c : aim.conformers)
		for (auto r : c)
		{
			auto gcrd = r.globalcrd;
			if (name == "ALL")
				for (auto it = gcrd.begin(); it != gcrd.end(); it++)
					crds.push_back(it->second);
			else if (name == "MC")
			{
				crds.push_back(gcrd["N"]);
				crds.push_back(gcrd["CA"]);
				crds.push_back(gcrd["C"]);
				crds.push_back(gcrd["O"]);
			}
			else
				crds.push_back(gcrd[name]);
		}
	return crds;
}

//
PYBIND11_MODULE(pb_proteinrep, m)
{
	py::class_<AAConformer>(m, "AAConformer")
	    .def(py::init<>())
	    .def("pb_make_aaconformer", overload_cast_<const std::vector<PdbRecord> &, std::vector<AAConformer> *, bool>()(&pb_make_aaconformer))
	    .def("pb_make_aaconformer", overload_cast_<NSPdesignseq::Rotamer *, const BackBoneSite *>()(&pb_make_aaconformer))
	    .def_readwrite_static("sidechainatoms", &AAConformer::sidechainatoms)
	    .def_readwrite_static("backboneatoms", &AAConformer::backboneatoms)
	    .def_readwrite_static("mainchainatoms", &AAConformer::mainchainatoms)
	    .def_readwrite("residuename", &AAConformer::residuename)
	    .def_readwrite("chainid_or", &AAConformer::chainid_or)
	    .def_readwrite("residueid_or", &AAConformer::residueid_or)
	    .def_readwrite("insertionid_or", &AAConformer::insertionid_or)
	    .def_readwrite("atomlist", &AAConformer::atomlist)
	    .def_readwrite("globalcrd", &AAConformer::globalcrd)
	    .def_readwrite("localcrd", &AAConformer::localcrd)
	    .def_readwrite("localframe", &AAConformer::localframe)
	    .def("calclocalframe", &AAConformer::calclocalframe)
	    .def("moveglobalcrd", &AAConformer::moveglobalcrd)
	    .def("calclocalcrd", &AAConformer::calclocalcrd)
	    .def("getlocalcrd", overload_cast_<>()(&AAConformer::getlocalcrd))
	    .def("getlocalcrd_const", overload_cast_<>()(&AAConformer::getlocalcrd, py::const_))
	    .def("calcglobalcrd", overload_cast_<>()(&AAConformer::calcglobalcrd))
	    .def("calcglobalcrd", overload_cast_<const std::map<std::string, NSPgeometry::XYZ> &>()(&AAConformer::calcglobalcrd))
	    .def("getglobalcrd", overload_cast_<>()(&AAConformer::getglobalcrd))
	    .def("getglobalcrd_const", overload_cast_<>()(&AAConformer::getglobalcrd, py::const_))
	    .def("mainchaincrdcomplete", &AAConformer::mainchaincrdcomplete)
	    .def("sidechaincrdcomplete", &AAConformer::sidechaincrdcomplete)
	    .def("crdcomplete", &AAConformer::crdcomplete)
	    .def("isnaturalaa", &AAConformer::isnaturalaa)
	    .def("isaa", &AAConformer::isaa)
	    .def("removeoxtorother", &AAConformer::removeoxtorother)
	    .def("removeatomsnotlisted", &AAConformer::removeatomsnotlisted)
	    .def("connectedto", &AAConformer::connectedto)
	    .def("connectedfrom", &AAConformer::connectedfrom)
	    .def("connected", &AAConformer::connected)
	    .def("getmainchaincrd", &AAConformer::getmainchaincrd)
	    .def("getbackbonecrd", &AAConformer::getbackbonecrd)
	    .def("phi", &AAConformer::phi)
	    .def("psi", &AAConformer::psi)
	    .def("omega_n", &AAConformer::omega_n)
	    .def("omega_c", &AAConformer::omega_c)
	    .def("dihedral", &AAConformer::dihedral)
	    .def("make_pdbrecords", &AAConformer::make_pdbrecords)
	    .def("make_backbonesite", &AAConformer::make_backbonesite)
	    .def("rotamerconformer", &AAConformer::rotamerconformer)
	    .def("sidechainRMSD2_local", &sidechainRMSD2_local)
	    .def("RMSD2_global", &RMSD2_global)
	    .def("backboneRMSD2", &backboneRMSD2)
//	    .def("distance", &distance)
	    .def("isbackbone", &isbackbone)
	    .def("ismainchain", &ismainchain);
	// AAConformer is all wrapped.

	py::class_<AAConformersInModel>(m, "AAConformersInModel")
	    .def(py::init<>())
	    .def_readwrite("pdbmodelid", &AAConformersInModel::pdbmodelid)
	    .def_readwrite("mappdbkeyint", &AAConformersInModel::mappdbkeyint)
	    .def_readwrite("conformers", &AAConformersInModel::conformers)
	    .def_readwrite("altconformers", &AAConformersInModel::altconformers)
//	    .def("getsequence", &AAConformersInModel::getsequence)
            .def("readpdbfile", &AAConformersInModel::readpdbfile)
            .def("getconformers", &AAConformersInModel::getconformers)
	    .def("getcrdbyname", &getcrdbyname)
	    .def("conformersinchain", &AAConformersInModel::conformersinchain);
	// AAConformersInModel is all wrapped.

	py::class_<PdbRecord> pdbrecord(m, "PdbRecord");
	pdbrecord.def(py::init<>())
	    .def(py::init<const std::string &>())
	    .def("init", &PdbRecord::init)
	    .def("make_pdbrecord", &mkrecord)
	    .def("_repr_", &PdbRecord::toString)
	    .def_readwrite("label",&PdbRecord::label)
	    .def_readwrite("atomname",&PdbRecord::atomname)
	    .def_readwrite("namesymbol",&PdbRecord::namesymbol)
	    .def_readwrite("namemodifier",&PdbRecord::namemodifier)
	    .def_readwrite("residuename",&PdbRecord::residuename)
	    .def_readwrite("atomid",&PdbRecord::atomid)
	    .def_readwrite("chainid",&PdbRecord::chainid)
	    .def_readwrite("residueid",&PdbRecord::residueid)
	    .def_readwrite("conformerid",&PdbRecord::conformerid)
	    .def_readwrite("insertionid",&PdbRecord::insertionid)
	    .def_readwrite("x",&PdbRecord::x)
	    .def_readwrite("y",&PdbRecord::y)
	    .def_readwrite("bfactor",&PdbRecord::bfactor)
	    .def_readwrite("occupation",&PdbRecord::occupation)
	    .def_property("elementname",get_pdr_elementname,set_pdr_elementname);
	py::enum_<PdbRecord::Field>(pdbrecord, "Field")
	    .value("LABEL", PdbRecord::Field::LABEL)
	    .value("ATOMID", PdbRecord::Field::ATOMID)
	    .value("NAMESYMBOL", PdbRecord::Field::NAMESYMBOL)
	    .value("NAMEMODIFIER", PdbRecord::Field::NAMEMODIFIER)
	    .value("CONFORMERID", PdbRecord::Field::CONFORMERID)
	    .value("RESIDUENAME", PdbRecord::Field::RESIDUENAME)
	    .value("CHAINID", PdbRecord::Field::CHAINID)
	    .value("RESIDUEID", PdbRecord::Field::RESIDUEID)
	    .value("INSERTIONID", PdbRecord::Field::INSERTIONID)
	    .value("X", PdbRecord::Field::X)
	    .value("Y", PdbRecord::Field::Y)
	    .value("Z", PdbRecord::Field::Z)
	    .value("OCCUPATION", PdbRecord::Field::OCCUPATION)
	    .value("BFACTOR", PdbRecord::Field::BFACTOR)
	    .value("SEGMENT", PdbRecord::Field::SEGMENT)
	    .value("ELEMENTNAME", PdbRecord::Field::ELEMENTNAME);
	// PdbRecord is all wrapped.

	py::class_<PdbReader>(m, "PdbReader")
	    .def(py::init<>())
	    .def("reskey", &PdbReader::reskey)
            .def("readpdb", overload_cast_<const std::string &>()(&PdbReader::readpdb))
            .def("readpdb", overload_cast_<std::vector<std::string> &>()(&PdbReader::readpdb))
	    .def("getaminoacidsequence", &PdbReader::getaminoacidsequence)
	    .def("chainids", &PdbReader::chainids)
	    .def("records", overload_cast_<>()(&PdbReader::records))
	    .def("records_const", overload_cast_<>()(&PdbReader::records, py::const_))
            .def("addRecord", &PdbReader::addRecord)
	    .def("mappdbkeyint", &PdbReader::mappdbkeyint);
	// PdbReader is all wrapped.

	py::class_<MapPdbKeyInt>(m, "MapPdbKeyInt")
	    .def(py::init<const typename PdbReader::RecordsType &, const std::map<char,std::map<int,std::vector<char>>> &, const std::vector<char> &>())
	    .def("pdbChainID", &MapPdbKeyInt::pdbChainID)
	    .def("pdbResKey", &MapPdbKeyInt::pdbResKey)
	    .def("pdbResID", &MapPdbKeyInt::pdbResID)
	    .def("chainNumber", &MapPdbKeyInt::chainNumber)
	    .def("posiNumber", &MapPdbKeyInt::posiNumber)
	    .def("mapchainidint", &MapPdbKeyInt::mapchainidint);
	// MapPdbKeyInt is all wrapped.

/*
	py::class_<IntAtomKey>(m, "IntAtomKey")
	    .def("genKey", &IntAtomKey::genKey)
	    .def("chainNumber", &IntAtomKey::chainNumber)
	    .def("posiNumber", &IntAtomKey::posiNumber)
	    .def("atomNmNumber", &IntAtomKey::atomNmNumber)
	    .def("resnameCode", &IntAtomKey::resnameCode)
	    .def("rotCodeNumber", &IntAtomKey::rotCodeNumber)
	    .def("atomName", &IntAtomKey::atomName)
	    .def("residueName", &IntAtomKey::residueName)
	    .def("positivelyCharged", overload_cast_<const std::string &>()(&IntAtomKey::positivelyCharged))
	    .def("positivelyCharged", overload_cast_<Key>()(&IntAtomKey::positivelyCharged))
	    .def("negativelyCharged", overload_cast_<const std::string &>()(&IntAtomKey::negativelyCharged))
	    .def("negativelyCharged", overload_cast_<Key>()(&IntAtomKey::negativelyCharged))
	    .def("isCharged", overload_cast_<const std::string &>()(&IntAtomKey::isCharged))
	    .def("isCharged", overload_cast_<Key>()(&IntAtomKey::isCharged))
	    .def("isPolar", overload_cast_<const std::string &>()(&IntAtomKey::isPolar))
	    .def("isPolar", overload_cast_<Key>()(&IntAtomKey::isPolar))
	    .def("isAromatic", overload_cast_<const std::string &>()(&IntAtomKey::isAromatic))
	    .def("isAromatic", overload_cast_<Key>()(&IntAtomKey::isAromatic))
	    .def("isMainChain", overload_cast_<const std::string &>()(&IntAtomKey::MainChain))
	    .def("isMainChain", overload_cast_<Key>()(&IntAtomKey::isMainChain))
	    .def("isHBDonor", overload_cast_<const std::string &>()(&IntAtomKey::isHBDonor))
	    .def("isHBDonor", overload_cast_<Key>()(&IntAtomKey::isHBDonor))
	    .def("isHBAcceptor", overload_cast_<const std::string &>()(&IntAtomKey::isHBAcceptor))
	    .def("isHBAcceptor", overload_cast_<Key>()(&IntAtomKey::isHBAcceptor));
	// IntAtomKey is all wrapped.
*/	
	//added by zl
	py::class_<AminoAcidSeq> aminoacidseq(m,"AminoAcidSeq");
    aminoacidseq.def("name2code",overload_cast_<const std::string &>()(&AminoAcidSeq::name2code))
        ;

    py::class_<BackBoneState> backbonestate(m,"BackBoneState");
    backbonestate.def_readwrite("phi",&BackBoneState::phi)
        .def_readwrite("psi",&BackBoneState::psi)
        .def_readwrite("omiga",&BackBoneState::omiga)
        .def_readwrite("sai",&BackBoneState::sai)
        .def_readwrite("phipsistate",&BackBoneState::phipsistate)
        .def_readwrite("SSState",&BackBoneState::SSState)
        .def_readwrite("SASAState",&BackBoneState::SASAState)
        ;

    py::class_<SideChainState> sidechainstate(m,"SideChainState");
    sidechainstate.def_readwrite("resiudetype",&SideChainState::resiudetype)
        .def_readwrite("rotamerstate",&SideChainState::rotamerstate)
        .def_readwrite("torsions",&SideChainState::torsions)
        ;

    py::class_<ResidueState> residuestate(m,"ResidueState");
    residuestate.def_readwrite("backbonestate",&ResidueState::backbonestate)
        .def_readwrite("sidechainstate",&ResidueState::sidechainstate)
        ;
	
	m.def("residuestates",[](std::string filename){
        std::ifstream is(filename);
        return residuestates(is);
    });
	//
}

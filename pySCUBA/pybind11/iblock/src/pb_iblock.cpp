#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>

#include "iblock/aasequtil.h"
#include "iblock/intrctmol.h"
#include "iblock/intrctblck.h"
#include "iblock/backboneff.h"
#include "iblock/blcktopo.h"
#include <iblock/molsystmpara.h>
#include <iblock/intrctparam.h>

#include "iblock/test/SCUBASketch.cpp"
#include "dstl/randomengine.h"

#include <fstream>
#include <iostream>
#include <string>

namespace py = pybind11;
namespace intr=NSPintrct;
using namespace NSPintrct;
namespace ibk = NSPintrct;

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;


typedef std::array<int, 3> AtomIdx3D; //blcktopo.h
typedef BondedIntrct<2> Bond;
typedef BondedIntrct<3> Angle;
typedef BondedIntrct<4> ImpDih;
typedef std::array<AtomIdx3D, 4> Torsion;

typedef NSPdstl::Idx2D SiteIdx; //molmodeler.h

typedef std::map<int, std::set<int>> SelectedBlks;


// struct USRIntrctMol: intr::IntrctMol
// {
//     void usrwritepdb(std::string filename){
//         std::ofstream ofs(filename);
//         writepdb(ofs);
//         ofs.close();
//     }
// };

struct USRIntrctMol: intr::IntrctMol
{
    void usrwritepdb(std::string filename){
        std::ofstream ofs(filename);
        writepdb(ofs);
        ofs.close();
    }
};

struct SCUBASketchMainPAR
{
  /* data */
  std::vector<char> sseq;
  std::vector<std::vector<double>> sscrd;
  std::vector<int> sslen;
  std::vector<char> finalss;
  std::vector<std::string>ss_info;
  std::vector<std::string>ssEH;
  std::vector<std::string>ssE;
  std::vector<std::string>ssH;
  std::vector<std::string>ssC;
  std::vector<char> direction;
  int linkloop;
  void readss(string file){
    ifstream infile;
    // std::cout<<"start file read"<<std::endl;
    infile.open(file.data());
    // std::cout<<"file is read"<<std::endl;
    assert(infile.is_open());
    string s;
    sseq.clear();
    sscrd.clear();
    // std::cout<<"getline"<<std::endl;
    while(getline(infile,s))
    {
        std::vector<double> temp_crd_start;
        std::vector<double> temp_crd_direc;
        std::vector<int> len_range;
        std::vector<std::string> words;
        NSPutils::split(s, words, ";");
        int flag=0;
        for (auto &word : words){
			std::vector<std::string> word1s;
			NSPutils::split(s, word1s, " ");
			for (auto &word1 : word1s){
				char c=word1[0];
				if ((c == 'H' or c == 'E') && flag==0 ){
					sseq.push_back(c);
					sseq.push_back(c);
				}
              else if ((flag==0)&& (c == 'N' or c == 'C')){
                direction.push_back(c);
              }
              else if (flag==0){
                len_range.push_back(stoi(word1));
              }
              if (flag==1){
                temp_crd_start.push_back(stod(word1));
              }
              if (flag ==2){
                temp_crd_direc.push_back(stod(word1));
              }
            }
            flag++;
        }

        if (!len_range.empty())
        {
            if (len_range.size() == 2)
            {
                int len = NSPdstl::RandomEngine<>::getinstance().intrng(len_range[0], len_range[1])();
                sslen.push_back(len);
            }
            else if (len_range.size() == 1)
            {
                int len = len_range[0];
                sslen.push_back(len);
            }
            else
            {
                std::cout << "Please input according to the rules" << std::endl;
                exit(1);
            }
        }
        sscrd.push_back(temp_crd_start);
        sscrd.push_back(temp_crd_direc);

    }
    infile.close();
  }
};
void writeSSinfo(std::string & OutputFile, std::vector<char>& finalss, std::vector<std::string>& ss_info, std::vector<std::string>& ssEH, std::vector<std::string>& ssE, std::vector<std::string>& ssH, std::vector<std::string>& ssC,int i) {
	std::string filename2 = OutputFile + std::to_string(i) + "_final_ss.txt";
	std::ofstream ofs1(filename2);
	for (auto c : finalss) ofs1 << c;
	ofs1 << std::endl;
	for (auto s : ss_info) ofs1 << s;
	ofs1 << std::endl;
	for (auto eh : ssEH) ofs1 << eh;
	ofs1 << std::endl;
	for (auto e : ssE) ofs1 << e;
	ofs1 << std::endl;
	for (auto h : ssH) ofs1 << h;
	ofs1 << std::endl;
	for (auto c : ssC) ofs1 << c;
	ofs1 << std::endl;
	ofs1.close();
}

void buildsketch(SketchPar& read_result) {
	std::vector<char> sseq;
	std::vector<std::vector<double>> sscrd;
	std::vector<int> sslen;
	std::vector<char> finalss;
	std::vector<std::string>ss_info;
	std::vector<std::string>ssEH;
	std::vector<std::string>ssE;
	std::vector<std::string>ssH;
	std::vector<std::string>ssC;
	std::vector<char> direction;

	int i = 0;
	while (i < read_result.Gennumber) {
		IntrctMol imol;
		MolModeler modeler;
		finalss.clear();
		ss_info.clear();
		ssEH.clear();
		ssE.clear();
		ssH.clear();
		ssC.clear();
		readss(read_result.SketchFile, sseq, sscrd, sslen, direction);
		if (i == 0) {
			NSPdstl::RandomEngine<>::getinstance().reseed(read_result.randomseed);
		}
		modeler.settarget(imol);
		modeler.buildssele(sseq, sscrd, sslen, read_result.linkloop, finalss, ss_info, ssEH, ssE, ssH, ssC, direction);
		std::string filename = read_result.OutputFile + std::to_string(i) + ".pdb";
		std::ofstream ofs(filename);
		imol.writepdb(ofs);
		ofs.close();
		std::string filename2 = read_result.OutputFile + std::to_string(i) + "_final_ss.txt";
		std::ofstream ofs1(filename2);
		for (auto c : finalss) ofs1 << c;
		ofs1 << std::endl;
		for (auto s : ss_info) ofs1 << s;
		ofs1 << std::endl;
		for (auto eh : ssEH) ofs1 << eh;
		ofs1 << std::endl;
		for (auto e : ssE) ofs1 << e;
		ofs1 << std::endl;
		for (auto h : ssH) ofs1 << h;
		ofs1 << std::endl;
		for (auto c : ssC) ofs1 << c;
		ofs1 << std::endl;
		ofs1.close();
		i++;
	}
}





PYBIND11_MODULE(pb_iblock, m) {
	//m.def("writeSSinfo", &writeSSinfo);
	m.def("buildsketch", &buildsketch);

//aasequtil.h
	m.def("getseq3l", &getseq3l);
	m.def("getseq1l", &getseq1l);

	//backboneff.h
	py::class_<NSPintrct::BackBoneFF>(m, "BackBoneFF")
		.def(py::init<>())
		.def_readwrite("natomspersites", &NSPintrct::BackBoneFF::natomspersites)
		.def_readwrite("b0_nca", &NSPintrct::BackBoneFF::b0_nca)
		.def_readwrite("kb_nca", &NSPintrct::BackBoneFF::kb_nca)
		.def_readwrite("b0_cac", &NSPintrct::BackBoneFF::b0_cac)
		.def_readwrite("kb_cac", &NSPintrct::BackBoneFF::kb_cac)
		.def_readwrite("b0_co", &NSPintrct::BackBoneFF::b0_co)
		.def_readwrite("kb_co", &NSPintrct::BackBoneFF::kb_co)
		.def_readwrite("b0_cn", &NSPintrct::BackBoneFF::b0_cn)
		.def_readwrite("kb_cn", &NSPintrct::BackBoneFF::kb_cn)
		.def_readwrite("t0_ncac", &NSPintrct::BackBoneFF::t0_ncac)
		.def_readwrite("kt_ncac", &NSPintrct::BackBoneFF::kt_ncac)
		.def_readwrite("t0_caco", &NSPintrct::BackBoneFF::t0_caco)
		.def_readwrite("kt_caco", &NSPintrct::BackBoneFF::kt_caco)
		.def_readwrite("t0_cacn", &NSPintrct::BackBoneFF::t0_cacn)
		.def_readwrite("kt_cacn", &NSPintrct::BackBoneFF::kt_cacn)
		.def_readwrite("t0_ocn", &NSPintrct::BackBoneFF::t0_ocn)
		.def_readwrite("kt_ocn", &NSPintrct::BackBoneFF::kt_ocn)
		.def_readwrite("t0_cnca", &NSPintrct::BackBoneFF::t0_cnca)
		.def_readwrite("kt_cnca", &NSPintrct::BackBoneFF::kt_cnca)
		.def_readwrite("p_cncao", &NSPintrct::BackBoneFF::p_cncao)
		.def_readwrite("kp_cncao", &NSPintrct::BackBoneFF::kp_cncao)
		.def_readwrite("p_cacnca", &NSPintrct::BackBoneFF::p_cacnca)
		.def_readwrite("kp_cacnca", &NSPintrct::BackBoneFF::kp_cacnca)
		.def_readwrite("p_ocnca", &NSPintrct::BackBoneFF::p_ocnca)
		.def_readwrite("kp_ocnca", &NSPintrct::BackBoneFF::kp_ocnca)
		.def_readwrite("sigma_n", &NSPintrct::BackBoneFF::sigma_n)
		.def_readwrite("sigma_nhb", &NSPintrct::BackBoneFF::sigma_nhb)
		.def_readwrite("sigma_ca", &NSPintrct::BackBoneFF::sigma_ca)
		.def_readwrite("sigma_c", &NSPintrct::BackBoneFF::sigma_c)
		.def_readwrite("sigma_o", &NSPintrct::BackBoneFF::sigma_o)
		.def_readwrite("sigma_ohb", &NSPintrct::BackBoneFF::sigma_ohb)
		.def_readwrite("eps", &NSPintrct::BackBoneFF::eps);

	//blcktopo.h

	//m.def("tostring", overload_cast_<const std::array<AtomIdx3D, NBODY> &>() (&NSPintrct::tostring));//模板全局变量报错,不用直接写模板
	//py::class_<NSPintrct::BondedIntrct>(m, "BondedIntrct")
		//.def(py::init<std::array<AtomIdx3D, NBODY>&, double&, double &>())
		//.def_readwrite("v0", &NSPintrct::BondedIntrct::v0)
		//.def_readwrite("k", &NSPintrct::BondedIntrct::k)
		//.def_readwrite("aidx", &NSPintrct::BondedIntrct::aidx)
		//.def("energy", &NSPintrct::BondedIntrct::energy);

	//m.def("Bond",&NSPintrct::Bond)
	//m.def("Angle", &NSPintrct::Angle)
	//m.def("ImpDih", &NSPintrct::ImpDih)
	//m.def("Torsion", &NSPintrct::Torsion)

	py::class_<NSPintrct::AtomIP>(m, "AtomIP")
		.def(py::init<int, double, double, double, bool, bool, bool, bool, bool>())
		.def("pairenergy", &NSPintrct::AtomIP::pairenergy)
		.def_readwrite("sigma", &NSPintrct::AtomIP::sigma)
		.def_readwrite("eps", &NSPintrct::AtomIP::eps)
		.def_readwrite("sigmahb", &NSPintrct::AtomIP::sigmahb)
		.def_readwrite("sterictype", &NSPintrct::AtomIP::sterictype)
		.def_readwrite("hbdonor", &NSPintrct::AtomIP::hbdonor)
		.def_readwrite("hbacceptor", &NSPintrct::AtomIP::hbacceptor)
		.def_readwrite("ismainchain", &NSPintrct::AtomIP::ismainchain)
		.def_readwrite("issidechain", &NSPintrct::AtomIP::issidechain)
		.def_readwrite("isprotein", &NSPintrct::AtomIP::hbacceptor)
		;

	py::class_<BlckTopo> blcktopo(m, "BlckTopo");
		//blcktopo.def_readwrite_static("CD_PRO", &NSPintrct::BlckTopo::CD_PRO)
	blcktopo.def_readwrite("blckname", &NSPintrct::BlckTopo::blckname)
		.def_static("getblcktopo_std", &NSPintrct::BlckTopo::getblcktopo_std)
		.def_static("getblcktopo_nonstd", &NSPintrct::BlckTopo::getblcktopo_nonstd)
		.def_readwrite("natoms", &NSPintrct::BlckTopo::natoms)
		.def_readwrite("atomnames", &NSPintrct::BlckTopo::atomnames)
		.def_readwrite("atomseqs", &NSPintrct::BlckTopo::atomseqs)
		.def_readwrite("bonds", &NSPintrct::BlckTopo::bonds)
		.def_readwrite("angles", &NSPintrct::BlckTopo::angles)
		.def_readwrite("impdihs", &NSPintrct::BlckTopo::impdihs)
		.def_readwrite("pepimpdih1", &NSPintrct::BlckTopo::pepimpdih1)
		.def_readwrite("pepimpdih2", &NSPintrct::BlckTopo::pepimpdih2)
		.def_readwrite("torsions", &NSPintrct::BlckTopo::torsions)
		.def_readwrite("atomips", &NSPintrct::BlckTopo::atomips)
		.def_readwrite("excld", &NSPintrct::BlckTopo::excld)
		.def_readwrite("nblist_internal_M", &NSPintrct::BlckTopo::nblist_internal_M)
		.def_readwrite("nblist_internal_S", &NSPintrct::BlckTopo::nblist_internal_S)
		.def("atomindex", &NSPintrct::BlckTopo::atomindex)
		.def("tostring", &NSPintrct::BlckTopo::tostring)
		.def("less", &NSPintrct::less);
	//m.def("setexcld", &setexcld);
	
	py::enum_<BlckTopo::BLCKPOSI>(blcktopo, "BLCKPOSI")
		.value("NTEMR", BlckTopo::BLCKPOSI::NTEMR)
		.value("CTEMR", BlckTopo::BLCKPOSI::CTEMR)
		.value("MIDDLE", BlckTopo::BLCKPOSI::MIDDLE);

		
	
	//blkselector.h
	py::class_<NSPintrct::BlkSelector>(m, "BlkSelector")
		.def(py::init<>())
		.def(py::init<const std::string&>())
		.def_readwrite("selectedblks", &NSPintrct::BlkSelector::selectedblks)
		.def("clear", &NSPintrct::BlkSelector::clear)
		.def("addselection", &NSPintrct::BlkSelector::addselection)
		.def("selectatoms", &NSPintrct::BlkSelector::selectatoms)
		.def("selectmcatoms", &NSPintrct::BlkSelector::selectmcatoms)
		.def("selectscatoms", &NSPintrct::BlkSelector::selectscatoms);

	//intrctmol.h
	m.def("addchain", &NSPintrct::IntrctMol::addchain);
	m.def("replaceblck", &NSPintrct::IntrctMol::replaceblck);
	m.def("pasteblck", &NSPintrct::IntrctMol::pasteblck);
	m.def("changecrds", overload_cast_<const std::vector<NSPgeometry::XYZ>&>()(&NSPintrct::IntrctMol::changecrds));
	m.def("changecrds", overload_cast_<const std::vector<NSPgeometry::XYZ>&>()(&NSPintrct::IntrctMol::changecrds, py::const_));
	m.def("nchains", &NSPintrct::IntrctMol::nchains);
	m.def("nresidues", &NSPintrct::IntrctMol::nresidues);
	m.def("natoms", overload_cast_<>()(&NSPintrct::IntrctMol::natoms, py::const_));
	m.def("natoms", overload_cast_<>()(&NSPintrct::IntrctMol::natoms));
	m.def("atomname", &NSPintrct::IntrctMol::atomname);
	m.def("residuename", &NSPintrct::IntrctMol::residuename);
	m.def("residueofatoms", overload_cast_<>()(&NSPintrct::IntrctMol::residueofatoms, py::const_));
	m.def("residueofatoms", overload_cast_<>()(&NSPintrct::IntrctMol::residueofatoms));
	m.def("atomindices", &NSPintrct::IntrctMol::atomindices);
	m.def("translatechaincenter", &NSPintrct::IntrctMol::translatechaincenter);
	m.def("recollectcrds_all", overload_cast_<>()(&NSPintrct::IntrctMol::recollectcrds_all, py::const_));
	m.def("recollectcrds_all", overload_cast_<>()(&NSPintrct::IntrctMol::recollectcrds_all));
	m.def("atomindices", &NSPintrct::IntrctMol::atomindices);
	m.def("getcrds_all_d", overload_cast_<>()(&NSPintrct::IntrctMol::getcrds_all_d, py::const_));
	m.def("getcrds_all_d", overload_cast_<>()(&NSPintrct::IntrctMol::getcrds_all_d));
	m.def("isfixedvec", overload_cast_<>()(&NSPintrct::IntrctMol::isfixedvec, py::const_));
	m.def("isfixedvec", overload_cast_<>()(&NSPintrct::IntrctMol::isfixedvec));
	m.def("changeactivemodeveryblck", &NSPintrct::IntrctMol::changeactivemodeveryblck);
	m.def("fixeveryblck", &NSPintrct::IntrctMol::fixeveryblck);
	m.def("activateeveryblck", &NSPintrct::IntrctMol::activateeveryblck);
	m.def("specifyactiveblks", &NSPintrct::IntrctMol::specifyactiveblks);
	m.def("calcphipsicodes", &NSPintrct::IntrctMol::calcphipsicodes);
	m.def("calcsscodes", &NSPintrct::IntrctMol::calcsscodes);
	m.def("sscodestrings", &NSPintrct::IntrctMol::specifyfixedblks);
	m.def("calcconformercodes", &NSPintrct::IntrctMol::calcconformercodes);
	m.def("clear_energies", &NSPintrct::IntrctMol::clear_energies);
	m.def("sum_energies", &NSPintrct::IntrctMol::sum_energies);
	m.def("sum_energies", [](NSPintrct::IntrctMol &mol){
		std::vector<NSPgeometry::XYZ> dedx;
		mol.my_intrctpara->phipsicodes_updateonly = false;				
		mol.my_intrctpara->sscodes_updateonly = false;				
		mol.my_intrctpara->calc_nblists = true;				
		mol.my_intrctpara->calc_sscodes = true;				
		// std::vector<NSPgeometry::XYZ> xyz = NSPgeometry::doublevtoXYZv(state_->crd);
		std::vector<NSPgeometry::XYZ> xyz = mol.recollectcrds_all();
		mol.forces_all(xyz, dedx);
		std::array<double,IntrctBlck::ENESIZE>  energies;
		double etot_=mol.sum_energies(energies);
		std::vector<double> totene;
		totene.push_back(etot_);
		for(auto e:energies) totene.push_back(e);
		return totene;
	});
	m.def("setintrctpara",[](NSPintrct::IntrctMol &mol, IntrctPara ipara){
		mol.my_intrctpara=std::shared_ptr<IntrctPara>(new IntrctPara(mkrescaledparam(ipara)));
	});

	m.def("mainchainwindowactive", &NSPintrct::IntrctMol::mainchainwindowactive);
	m.def("forces_all", overload_cast_<const IntrctPara&,std::vector<NSPgeometry::XYZ>&>()(&NSPintrct::IntrctMol::forces_all, py::const_));
	m.def("forces_all", overload_cast_<std::vector<NSPgeometry::XYZ>&>()(&NSPintrct::IntrctMol::forces_all, py::const_));
	m.def("forces_all", overload_cast_<const std::vector<NSPgeometry::XYZ>&,std::vector<NSPgeometry::XYZ>&>()(&NSPintrct::IntrctMol::forces_all, py::const_));
	m.def("forces_all", overload_cast_<const IntrctPara&,const std::vector<NSPgeometry::XYZ>&,std::vector<NSPgeometry::XYZ>&>()(&NSPintrct::IntrctMol::forces_all, py::const_));
	m.def("forces_cov", &NSPintrct::IntrctMol::forces_cov);
	m.def("mksteric_nblists", &NSPintrct::IntrctMol::mksteric_nblists);
	m.def("forces_steric", &NSPintrct::IntrctMol::forces_steric);
	m.def("forces_phipsi", &NSPintrct::IntrctMol::forces_phipsi);
	m.def("forces_localstr", &NSPintrct::IntrctMol::forces_localstr);
	m.def("forces_sitepair", &NSPintrct::IntrctMol::forces_sitepair);
	m.def("forces_rotamer", &NSPintrct::IntrctMol::forces_rotamer);
	m.def("forces_scpacking", &NSPintrct::IntrctMol::forces_scpacking);
	m.def("forces_localhb", &NSPintrct::IntrctMol::forces_localhb);
	m.def("ssregions", &NSPintrct::IntrctMol::ssregions);
	m.def("writepdb", &NSPintrct::IntrctMol::writepdb);
	m.def("checkblckindices", &NSPintrct::IntrctMol::checkblckindices);
	m.def("checkaoffsets", &NSPintrct::IntrctMol::checkaoffsets);
	m.def("checkbsinchains", &NSPintrct::IntrctMol::checkbsinchains);
	m.def("checkatomsfixed", &NSPintrct::IntrctMol::checkatomsfixed);
	m.def("checkscinchains", &NSPintrct::IntrctMol::checkscinchains);
	m.def("checkcrds", &NSPintrct::IntrctMol::checkcrds);
	m.def("sum_energies", &NSPintrct::IntrctMol::sum_energies);
	m.def("make_molsystm", &NSPintrct::make_molsystm);
	

	py::class_<NSPintrct::IntrctMol>(m, "IntrctMol")
		.def(py::init<>())
		.def(py::init<const NSPintrct::IntrctMol&>())
		.def(py::init<const std::vector<std::string>>())
		.def(py::init<const std::vector<std::vector<std::string>>>())
		.def(py::init<const std::vector<std::vector<NSPproteinrep::AAConformer>>&>())
		.def(py::init<const NSPproteinrep::AAConformersInModel&>())
		.def("setup", overload_cast_<const std::vector<std::vector<std::string>>&>()(&NSPintrct::IntrctMol::setup))
		.def("setup", overload_cast_<const std::vector<std::vector<NSPproteinrep::AAConformer>>&, bool >()(&NSPintrct::IntrctMol::setup))
		.def("setup", overload_cast_<const std::vector<std::string>&>()(&NSPintrct::IntrctMol::setup))
		.def("getmolsystm", overload_cast_<>()(&NSPintrct::IntrctMol::getmolsystm))
		.def("getmolsystm", overload_cast_<>()(&NSPintrct::IntrctMol::getmolsystm, py::const_))
		.def("getblck", overload_cast_<const NSPdstl::Idx2D&>()(&NSPintrct::IntrctMol::getblck))
		.def("getblck", overload_cast_<const NSPdstl::Idx2D&>()(&NSPintrct::IntrctMol::getblck, py::const_))
		.def("__call__", overload_cast_<const NSPdstl::Idx2D&>()(&NSPintrct::IntrctMol::operator()))
		.def("__call__", overload_cast_<const NSPdstl::Idx2D&>()(&NSPintrct::IntrctMol::operator(), py::const_))
		.def("mappdbkeyint", overload_cast_<>()(&NSPintrct::IntrctMol::mappdbkeyint))
		.def("mappdbkeyint", overload_cast_<>()(&NSPintrct::IntrctMol::mappdbkeyint, py::const_))

		//added by zl
		.def("writepdb",[](IntrctMol &mol, std::string filename){
            std::ofstream ofs(filename);
            mol.writepdb(ofs);
        })
		.def("natoms", overload_cast_<>()(&NSPintrct::IntrctMol::natoms))
		.def("specifyactiveblks",[](IntrctMol &mol, std::vector<std::pair<int,std::set<int>>> mcvec, std::vector<std::pair<int,std::set<int>>> scvec){
			BlkSelector::SelectedBlks mcsel,scsel; //std::map<int,std::set<int>>
			for(int n=0;n<mcvec.size();n++){
				mcsel[mcvec[n].first]=mcvec[n].second;
			}
			for(int n=0;n<scvec.size();n++){
				scsel[scvec[n].first]=scvec[n].second;
			}
			mol.specifyactiveblks(mcsel,scsel);
		})
		//added by zl

		;

	py::class_<NSPintrct::IntrctBlck>(m, "IntrctBlck")
		.def(py::init<>())
	;

	//molmodeler.h
	//m.def("mkimol_newloop", &NSPintrct::MolModeler::mkimol_newloop); //copyclass
	//m.def("mergechains", &NSPintrct::MolModeler::mergechains);    //?与上面相同但是不能编译
	//m.def("splitchain", &NSPintrct::MolModeler::splitchain);
	m.def("ssregions", &NSPintrct::ssregions);

	py::class_<NSPintrct::MolModeler>(m, "MolModeler")
        .def(py::init<>())
		.def("settarget", &NSPintrct::MolModeler::settarget)
		.def("getbackbonesites", overload_cast_<>()(&NSPintrct::MolModeler::getbackbonesites))
		.def("getbackbonesites", overload_cast_<>()(&NSPintrct::MolModeler::getbackbonesites, py::const_))
		.def("checkclash", &NSPintrct::MolModeler::checkclash)
		.def("setneigborcutoff", &NSPintrct::MolModeler::setneigborcutoff)
		.def("changeresidue", &NSPintrct::MolModeler::changeresidue)
		.def("mkimol_newloop", &NSPintrct::MolModeler::mkimol_newloop) //copyclass
		//.def("mergechains", &NSPintrct::MolModeler::mergechains)//?与上面相同但是不能编译
		.def("splitchain", &NSPintrct::MolModeler::splitchain)
		.def("buildnewchain", &NSPintrct::MolModeler::buildnewchain)
		.def("buildssele", &NSPintrct::MolModeler::buildssele)
		.def("movechaincenter", &NSPintrct::MolModeler::movechaincenter)
		.def("replaceloop", &NSPintrct::MolModeler::replaceloop)
		.def("outputloop", &NSPintrct::MolModeler::outputloop)
		.def("newloop", &NSPintrct::MolModeler::newloop)
		.def("getnearbysites", overload_cast_<const SiteIdx&>()(&NSPintrct::MolModeler::getnearbysites))
		.def("getnearbysites", overload_cast_<const SiteIdx&>()(&NSPintrct::MolModeler::getnearbysites, py::const_))
		.def("getphipsi", overload_cast_<const SiteIdx&>()(&NSPintrct::MolModeler::getphipsi))
		.def("getphipsi", overload_cast_<const SiteIdx&>()(&NSPintrct::MolModeler::getphipsi, py::const_))
		.def("getlocalframe", overload_cast_<const SiteIdx&>()(&NSPintrct::MolModeler::getlocalframe))
		.def("getlocalframe", overload_cast_<const SiteIdx&>()(&NSPintrct::MolModeler::getlocalframe, py::const_))
		.def("nclashes", &NSPintrct::MolModeler::nclashes)
		.def("buildsidechain", overload_cast_<const SiteIdx&,const std::string&>()(&NSPintrct::MolModeler::buildsidechain))
		.def("buildsidechain", overload_cast_<const SiteIdx&,const std::string&, bool>()(&NSPintrct::MolModeler::buildsidechain))
		// .def("mkfloatingblck", overload_cast_<const SiteIdx&, NSPdesignseq::Rotamer*>()(&NSPintrct::MolModeler::mkfloatingblck))
		// .def("mkfloatingblck", overload_cast_<const SiteIdx&,const std::string, const NSPproteinrep::BackBoneSite&>()(&NSPintrct::MolModeler::mkfloatingblck, py::const_))
		.def("getloopconfmaker", &NSPintrct::MolModeler::getloopconfmaker)
		.def("mkloopreservoir", &NSPintrct::MolModeler::mkloopreservoir)

		//added by zl
		.def("settarget",[](MolModeler &modeler,IntrctMol &target){
            modeler.settarget(target);
        })
        .def("changeresidue",[](MolModeler &modeler,std::pair<int,int> chainid_resid, std::string resname){
            modeler.changeresidue(MolModeler::SiteIdx{chainid_resid.first,chainid_resid.second},resname);
        })
		.def("mkfloatingblck",[](MolModeler &modeler,std::pair<int,int> site,std::string restype, NSPproteinrep::BackBoneSite &bs){
            return modeler.mkfloatingblck(MolModeler::SiteIdx{site.first,site.second},restype,bs);
        })
		//added by zl
		;


	//zjw
	//molsystmpara.h
	py::class_<ibk::MolSystmPara>(m, "MolSystmPara")
		.def(py::init<>())
		.def(py::init<const std::vector<std::string>&>())
		.def_readwrite("sequence1l", &ibk::MolSystmPara::sequence1l)
		.def_readwrite("pdbstart", &ibk::MolSystmPara::pdbstart)
		.def_readwrite("jobname", &ibk::MolSystmPara::jobname)
		.def_readwrite("fixedmainchainresidues", &ibk::MolSystmPara::fixedmainchainresidues)
		.def_readwrite("fixedresidues", &ibk::MolSystmPara::fixedresidues)
		.def_readwrite("activeresidues", &ibk::MolSystmPara::activeresidues)
		.def_readwrite("sidechainactiveresidues", &ibk::MolSystmPara::sidechainactiveresidues)
		.def_readwrite("softsidechainresidues", &ibk::MolSystmPara::softsidechainresidues);
	m.def("makemolsystmparam", &ibk::makemolsystmparam);
	//intrctparam.h 
	py::class_<ibk::Results>(m, "Results")
		.def(py::init<>())
		.def_static("ofstreams", &ibk::Results::ofstreams);
	py::class_<ibk::IntrctPara>(m, "IntrctPara")
		.def(py::init<>())
		.def(py::init<const std::vector<std::string>&>())
		.def_readwrite("enedetails", &ibk::IntrctPara::enedetails)
		.def_readwrite("phipsicodes_updateonly", &ibk::IntrctPara::phipsicodes_updateonly)
		.def_readwrite("sscodes_updateonly", &ibk::IntrctPara::sscodes_updateonly)
		.def_readwrite("calc_sscodes", &ibk::IntrctPara::calc_sscodes)
		.def_readwrite("calc_nblists", &ibk::IntrctPara::calc_nblists)
		.def_readwrite("jobname", &ibk::IntrctPara::jobname)
		.def_readwrite("mm_neighbor_cutoff2", &ibk::IntrctPara::mm_neighbor_cutoff2)
		.def_readwrite("ms_neighbor_cutoff2", &ibk::IntrctPara::ms_neighbor_cutoff2)
		.def_readwrite("ss_neighbor_cutoff2", &ibk::IntrctPara::ss_neighbor_cutoff2)
		.def_readwrite("pl_neighbor_cutoff2", &ibk::IntrctPara::pl_neighbor_cutoff2)
		.def_readwrite("ll_neighbor_cutoff2", &ibk::IntrctPara::ll_neighbor_cutoff2)
		.def_readwrite("max_cutoff2", &ibk::IntrctPara::max_cutoff2)
		.def_readwrite("weight_coval", &ibk::IntrctPara::weight_coval)
		.def_readwrite("weight_steric", &ibk::IntrctPara::weight_steric)
		.def_readwrite("weight_phipsi", &ibk::IntrctPara::weight_phipsi)
		.def_readwrite("weight_localstr", &ibk::IntrctPara::weight_localstr)
		.def_readwrite("weight_sitepair", &ibk::IntrctPara::weight_sitepair)
		.def_readwrite("weight_rotamer", &ibk::IntrctPara::weight_rotamer)
		.def_readwrite("weight_scpacking", &ibk::IntrctPara::weight_scpacking)
		.def_readwrite("weight_localhb", &ibk::IntrctPara::weight_localhb);
	m.def("makeintrctparam", &ibk::makeintrctparam);
	m.def("mkrescaledparam", &ibk::mkrescaledparam);


    py::class_<USRIntrctMol,intr::IntrctMol>(m, "USRIntrctMol")
        .def(py::init<>())
        .def("usrwritepdb", &USRIntrctMol::usrwritepdb);

    py::class_<SketchPar>(m, "SketchPar")
        .def(py::init<>())
        .def(py::init<const std::vector<std::string> & >())
        .def_readwrite("SketchFile",&SketchPar::SketchFile)
        .def_readwrite("linkloop",&SketchPar::linkloop)
        .def_readwrite("OutputFile",&SketchPar::OutputFile)
        .def_readwrite("Gennumber",&SketchPar::Gennumber)
        .def_readwrite("randomseed",&SketchPar::randomseed)
        .def("readss", &readss)
        .def("readparameters", &readparameters);
    
    py::class_<SCUBASketchMainPAR>(m, "SCUBASketchMainPAR")
        .def(py::init<>())
        .def_readwrite("sseq",&SCUBASketchMainPAR::sseq)
        .def_readwrite("sscrd",&SCUBASketchMainPAR::sscrd)
        .def_readwrite("sslen",&SCUBASketchMainPAR::sslen)
        .def_readwrite("finalss",&SCUBASketchMainPAR::finalss)
        .def_readwrite("ss_info",&SCUBASketchMainPAR::ss_info)
        .def_readwrite("ssEH",&SCUBASketchMainPAR::ssEH)
        .def_readwrite("ssE",&SCUBASketchMainPAR::ssE)
        .def_readwrite("ssH",&SCUBASketchMainPAR::ssH)
        .def_readwrite("ssC",&SCUBASketchMainPAR::ssC)
        .def_readwrite("direction",&SCUBASketchMainPAR::direction)
        .def("readss",&SCUBASketchMainPAR::readss);
    
    
};
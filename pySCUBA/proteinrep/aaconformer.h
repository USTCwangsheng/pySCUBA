/*
 *aaconformer.h
 *
 *  Created on: 2017年6月23日
 *      Author: hyliu
 */

#ifndef PROTEINREP_AACONFORMER_H_
#define PROTEINREP_AACONFORMER_H_
#include "geometry/xyz.h"
#include "geometry/localframe.h"
#include "proteinrep/pdbreader.h"
#include "designseq/ProteinRep.h"
#include <memory>

namespace NSPproteinrep {
struct BackBoneSite;
struct AAConformer {
	static std::map<std::string, std::vector<std::string>> sidechainatoms;
	static std::vector<std::string> backboneatoms;
	static std::vector<std::string> mainchainatoms;
	std::string residuename { "" };
	char chainid_or{'A'}; ///<original chainid id if from a PDB file
	int residueid_or{0}; ///<original residue id if from a PDB file
	char insertionid_or{' '};///original inseretionid from a PDB file
	//list of atoms, preserves the order of atoms in the PDB records
	std::vector<std::string> atomlist;
	std::map<std::string, NSPgeometry::XYZ> globalcrd;
	std::map<std::string, NSPgeometry::XYZ> localcrd;
	NSPgeometry::LocalFrame localframe;

	template<class Archive>
	    void serialize(Archive & ar, const unsigned int version)
	    {
	        ar & residuename;
	        ar & chainid_or;
	        ar & residueid_or;
	        ar & atomlist;
	        ar & globalcrd;
	        ar & localcrd;
	    }
	static bool ismainchain(const std::string &atomname){
		for(auto &s:mainchainatoms) if(s==atomname) return true;
		return false;
	}
	static NSPgeometry::LocalFrame calclocalframe(
			const std::map<std::string, NSPgeometry::XYZ> & gcrd,
			std::vector<std::string> frameatoms=std::vector<std::string>());
	std::map<std::string, NSPgeometry::XYZ> moveglobalcrd(
			const std::map<std::string,NSPgeometry::XYZ> &refgcrd,
			const std::vector<std::string> &frmatoms=std::vector<std::string>());
	std::map<std::string, NSPgeometry::XYZ> & calclocalcrd(
			const std::vector<std::string> &frameatoms=std::vector<std::string>());
	const std::map<std::string, NSPgeometry::XYZ> &getlocalcrd() const {
		return localcrd;
	}
	std::map<std::string, NSPgeometry::XYZ> &getlocalcrd() {
		return localcrd;
	}
	std::map<std::string, NSPgeometry::XYZ> & calcglobalcrd(
			const std::map<std::string, NSPgeometry::XYZ> &refglobal) {
		localframe = calclocalframe(refglobal);
		return calcglobalcrd();
	}
	std::map<std::string, NSPgeometry::XYZ> & calcglobalcrd();
	const std::map<std::string, NSPgeometry::XYZ> & getglobalcrd() const {
		return globalcrd;
	}
	std::map<std::string, NSPgeometry::XYZ> & getglobalcrd() {
		return globalcrd;
	}
	bool mainchaincrdcomplete()const {
		for (auto &m : mainchainatoms) {
			if (globalcrd.find(m) == globalcrd.end())
				return false;
		}
		return true;
	}
	bool sidechaincrdcomplete() const {
		if(!isnaturalaa()) return true;
		for (auto & a : sidechainatoms.at(residuename)) {
			if (globalcrd.find(a) == globalcrd.end())
				return false;
		}
		return true;
	}
	bool crdcomplete() const {
		if(!isnaturalaa()) return true;
		return mainchaincrdcomplete() && sidechaincrdcomplete();
	}
	bool isnaturalaa() const {
		return (sidechainatoms.find(residuename) != sidechainatoms.end());
	}
	bool isaa() const {
		return isnaturalaa() || mainchaincrdcomplete();
		}
	void removeoxtorother(const std::string & name="OXT"){
		globalcrd.erase(name);
		localcrd.erase(name);
	}
	int removeatomsnotlisted();
	bool connectedto(const AAConformer &cn) const {
		return ((globalcrd.at("C")-cn.globalcrd.at("N")).squarednorm()<=3.24);
	}
	bool connectedfrom(const AAConformer &cp) const{
		return cp.connectedto(*this);
	}
	bool connected(const AAConformer &c) const {
		return connectedto(c) || connectedfrom(c);
	}
	std::vector<NSPgeometry::XYZ> getmainchaincrd() const{
		std::vector<NSPgeometry::XYZ> res;
		for(auto &a:mainchainatoms) res.push_back(globalcrd.at(a));
		return res;
	}
	std::vector<NSPgeometry::XYZ> getbackbonecrd() const{
			std::vector<NSPgeometry::XYZ> res;
			for(auto &a:backboneatoms) res.push_back(globalcrd.at(a));
			return res;
		}
	double phi(const AAConformer &cp) const {
		return NSPgeometry::torsion(cp.globalcrd.at("C"), globalcrd.at("N"),
				globalcrd.at("CA"), globalcrd.at("C"));
	}
	double psi(const AAConformer &cn) const {
		return NSPgeometry::torsion(globalcrd.at("N"), globalcrd.at("CA"),
				globalcrd.at("C"), cn.globalcrd.at("N"));
	}
	double omega_n(const AAConformer &cp) const {
		return NSPgeometry::torsion(cp.globalcrd.at("CA"), cp.globalcrd.at("C"),
				globalcrd.at("N"), globalcrd.at("CA"));
	}
	double omega_c(const AAConformer &cn) const {
		return NSPgeometry::torsion(globalcrd.at("CA"), globalcrd.at("C"),
				cn.globalcrd.at("N"), cn.globalcrd.at("CA"));
	}
	double dihedral(const std::string &a1, const std::string &a2,
			const std::string &a3, const std::string &a4) {
		if (!globalcrd.empty()) {
			return NSPgeometry::torsion(globalcrd.at(a1), globalcrd.at(a2),
					globalcrd.at(a3), globalcrd.at(a4));
		} else {
			return NSPgeometry::torsion(localcrd.at(a1), localcrd.at(a2),
					localcrd.at(a3), localcrd.at(a4));
		}
	}
	std::vector<PdbRecord> make_pdbrecords(char chainid, int resid,
			char insertionid, int atomid0) const;
	BackBoneSite make_backbonesite() const;
	/**
	 * make an AAConformer in the closest rotamer state of the given conformer
	 */
	AAConformer rotamerconformer() const;
};

double sidechainRMSD2_local(const AAConformer &c1, const AAConformer &c2);
double RMSD2_global(const AAConformer &c1, const AAConformer &c2);
double backboneRMSD2(const AAConformer &c1, const AAConformer &c2);

double distance(const AAConformer &c1, const AAConformer &c2,
		const std::string & atom1 = "CA", const std::string &atom2 = "CA");
bool isbackbone(const std::string &atomname);
bool ismainchain(const std::string &atomname);
AAConformer make_aaconformer(const std::vector<PdbRecord> & residuerecords,
		std::vector<AAConformer> *altconformers = nullptr,bool replaceMSE=true);

/**
 * make an AAConformer in a given rotamer state
 * If site !=nullptr, the backbone atom coordinates in site will be used.
 */
AAConformer make_aaconformer(NSPdesignseq::Rotamer *rt,
		const BackBoneSite *site=nullptr);
struct AAConformersInModel {
	PDBModelID pdbmodelid;
	std::shared_ptr<MapPdbKeyInt> mappdbkeyint { nullptr };
	std::vector<std::vector<AAConformer>> conformers;
	std::vector<std::vector<std::vector<AAConformer>>>altconformers;
	std::string getsequence(int chainnumber) const;
	void readpdbfile(const std::string &filename);
	void getconformers(const PdbReader &pdb);
	const std::vector<AAConformer> & conformersinchain(char chainid=' ') const
	{
		if( chainid != ' ' && mappdbkeyint) {
			return conformers.at(mappdbkeyint->chainNumber(chainid));
		}
		else return conformers[0];
	};
};

}

#endif /* PROTEINREP_AACONFORMER_H_ */

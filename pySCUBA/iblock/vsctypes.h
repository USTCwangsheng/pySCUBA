/*
 * sidechainff.h
 *
 *  Created on: 2018年1月31日
 *      Author: hyliu
 */

#ifndef IBLOCK_SIDECHAINFF_H_
#define IBLOCK_SIDECHAINFF_H_
#include "iblock/enefunc1d.h"
#include "iblock/phipsicodes.h"
#include <vector>
#include <string>
#include <map>
#include <set>
#include <memory>
#include <cassert>
#define A2NM 0.1
//#define MAXEXCL  15
#define KBT 2.4942
#define KANG_FAC 3282.72
namespace NSPintrct{
/**
 * different ways of calculating sidechain packing (normal, soft...)
 */
class PackingEneFuncs{
public:
	EneFunc1D & getenefunc(int atype,int btype) const {
		return *(functable_[eidx(atype,btype)]);
	}
	void setup(int ntypes){
		ntypes_=ntypes;
		functable_.assign(ntypes*ntypes,nullptr);
	}
	void addfunction(int a,int b,std::shared_ptr<EneFunc1D> func){
		assert(a<ntypes_ && b<ntypes_);
		functable_[eidx(a,b)]=func;
	}
private:
	std::vector<std::shared_ptr<EneFunc1D>> functable_;
	int ntypes_{0};
	int eidx(int a, int b) const {return a*ntypes_+b;}
};
/**
 * atom type for sidechain packing energy calculation
 */
struct PackingAtomType{
	double radius;
	int hbtype;
	int aromatic;
	int nconnections;
	static bool hbond(int hb1,int hb2){
		if(hb1 == 1 ){
			if(hb2==2 || hb2==3) return true;
		} else if(hb2==1){
			if(hb1==2 || hb1==3) return true;
		}
		return false;
	}
};
/**
 * defines topology-like structures of standard residue types, which can be read from input files
 */
struct VSCType{
	/**
	 * returns a reference to the object defined for a given residue type
	 *
	 * On first call, will read an input data file and construct the objects based on the input data.
	 */
	static const VSCType & getVSCType(const std::string & resname, const std::string &filename=std::string());
	static std::map<std::string,VSCType> readVSCTypes(const std::string & filename);
	static std::vector<PackingAtomType>packingatomtypes;
	static PackingEneFuncs packingenefuncs;
	static PackingEneFuncs softpackingenefuncs;
	static std::map<std::string,int> stericatomtypes;
	static std::map<char,std::string> resnamefrom1letter;
	static std::map<std::string,std::set<std::string>> rotatablescatoms;
	static bool isrotatablescatom(const std::string & resname, const std::string & atomname) {
	    return rotatablescatoms.find(resname) != rotatablescatoms.end()
	            && rotatablescatoms.at(resname).find(atomname) != rotatablescatoms.at(resname).end();
	}
	static std::string resnameof( char letter){
		if(resnamefrom1letter.empty()) getVSCType("ALA");
		return resnamefrom1letter.at(letter);
	}
	static int getstericatomtype(const std::string &resname,const std::string &atomname){
		std::string fullname=resname+":"+atomname;
		auto it=stericatomtypes.find(fullname);
		if(it== stericatomtypes.end()){
			fullname=std::string("ANY:")+atomname;
			it=stericatomtypes.find(fullname);
		}
		if(it==stericatomtypes.end()) return -1;
		return  it->second;
	}
	int nscatoms; //*< number of sidechain atoms
	std::string resname;
	std::string pdbname;
	char oneletter{'\0'};
	std::vector<std::string> atomnames; //*<ordered side chain atom names
	std::vector<std::vector<std::pair<int,double>>> internalcrds; //*<internal coordinates sequentially defining the sidechain atom coordinates
	std::vector<int> rotameratoms; //*< atoms whose torsional angles defines the rotamer kai angles
	std::vector<std::pair<int,int>> newbonds;
	std::vector<double> b0;
	std::vector<double> kb0;
	std::vector<std::vector<int>> newangles;
	std::vector<double> a0;
	std::vector<double> ka0;
	std::vector<std::vector<int>> newimpdihs;
	std::vector<double> imp0;
	std::vector<double> kimp0;
	std::vector<std::vector<int>> newtorsions;
};

/**
 * sidechain atoms in a molecule system
 */
struct SCInChain {
	std::string restype;
	int poffset{-1}; //*<absolute atom index of the CB atom;
	int nscatoms{0};
	bool softpacking{false};
	std::vector<std::vector<int>> kaiatoms; //*<absolute atom index of atoms for calculation the kai angles
	SCInChain(){;}
	SCInChain(const std::string & name,
			const std::vector<std::vector<int>> &kais,int caposi,int nasc):restype(name),
					kaiatoms(kais),poffset(caposi),nscatoms(nasc){;}
};
std::vector<std::vector<SCInChain>> makescinchains(
		const IntrctMol &mol);
struct ConformerCode {
	PhiPsiCodes *phipsicodes;
	std::vector<double> sidechaintorsions;
	std::vector<std::vector<double>> torsioncodes;
	std::vector<std::vector<std::vector<DvDxi>>> dtorsioncodesdx;
	std::string restype;
	static std::vector<double> gettorsioncodes(double ang,
			const std::vector<DvDxi> &dadx, std::vector<std::vector<DvDxi>> *dcdx);
};
std::vector<ConformerCode> makeconformercodes(const std::vector<NSPgeometry::XYZ> &crds,
		const std::vector<SCInChain> &scinchains,
		std::vector<PhiPsiCodes> &phipsicodes);
double mcscpackingenergy(const std::vector<NSPgeometry::XYZ> &crds, const std::vector<int> &atomtypes,
		const BSInChain &mc,const SCInChain &sc,
		int sep,std::vector<DvDxi> *dedx,bool terminal=false);
double scscpackingenergy(const std::vector<NSPgeometry::XYZ> &crds,
		const std::vector<int> &atomtypes,
		const SCInChain &sc1, const SCInChain &sc2,
		int sep,std::vector<DvDxi> *dedx);
}



#endif /* IBLOCK_SIDECHAINFF_H_ */

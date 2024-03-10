/*
 * phipsicodes.h
 *
 *  Created on: 2019年12月4日
 *      Author: hyiu
 */

#ifndef IBLOCK_PHIPSICODES_H_
#define IBLOCK_PHIPSICODES_H_
#include "geometry/xyz.h"
#include "backbone/backbonesite.h"
#include <vector>
#include <cassert>

namespace NSPintrct{
/**
 * Backbone atoms in a molecule system
 */
struct BSInChain{
	BSInChain(){;}
	BSInChain(int n,int ca, int c, int o,bool a): nid(n),caid(ca),cid(c),oid(o),active(a){;}
	int nid{0},caid{1},cid{2},oid{3}; //*< absolute atom indices of backbone atoms
	int chainid{0}; //*< chainid of the  backbone position
	int resid{0}; //*< residue id (in chain) of the position
	bool active{true};
	std::vector<int> atomids() const{
		std::vector<int> res;
		res.push_back(nid);res.push_back(caid);res.push_back(cid);res.push_back(oid);
		return res;
	}
};
/**
 * determines whether any backbone atoms is active in a window around a central position
 */
bool windowactive(const std::vector<BSInChain> &bsinchain, int center,int halfw);

/**
 * binds atom ID with a derivative with respect to the atom coordinates
 */
typedef std::pair<int, NSPgeometry::XYZ> DvDxi;
inline DvDxi operator*(double c, const DvDxi &dvdxi){
	return std::make_pair(dvdxi.first,c*(dvdxi.second));
}
inline DvDxi operator+(const DvDxi & d1, const DvDxi &d2){
	assert(d1.first==d2.first);
	return std::make_pair(d1.first, d1.second+d2.second);
}
inline NSPgeometry::XYZ getxyz(const std::vector<double> &crd, int posi){
	int idx=3*posi;
	return NSPgeometry::XYZ(crd[idx],crd[idx+1],crd[idx+2]);
}

/**
 * phipsi and their encoded values (for input to a neuronetwork)
 */
struct PhiPsiCodes{
    static std::vector<double> gettriangcodes(double ang,const std::vector<DvDxi> &dadx,
    		std::vector<std::vector<DvDxi>>  *dcdx);
    PhiPsiCodes(){;}
    PhiPsiCodes(const std::vector<NSPgeometry::XYZ> &crd, std::vector<BSInChain>::const_iterator bs,
    		bool Nterm=false,bool Cterm=false);
    double phi{360.0};
    double psi{360.0};
    std::vector<DvDxi> dphidx; //*< derivatives of phi w.r.t. coordinates
    std::vector<DvDxi> dpsidx; //*< derivatives of psi w.r.t. coordinates
	std::vector<double> phicodes; //*<encoding vector for phi
	std::vector<double> psicodes; //*encoding vector for psi
	std::vector<std::vector<DvDxi>> dphicodesdx; //*< derivative of the encoded values, dimension: code_length X 4
	std::vector<std::vector<DvDxi>> dpsicodesdx; //*<derivatives of the encoded values
};

/**
* encodes secondary structures
*/
struct SSCode{
	std::vector<double> p3; //*<estimated probability for eacth SS type
	std::vector<std::vector<DvDxi>> dp3dx; //*< derivatives of the probabilities w.r.t atomic coordinates
	int ssid; //*< SS type ID
	char charcode()const {
		static std::string types{"HECT"}; //helix,strand,coil,terminus
		return types[ssid];
	}
};
class IntrctMol;
/**
 * extract backbone sites in the Molecular system
 */
std::vector<std::vector<BSInChain>> makebsinchains(const IntrctMol &mol);

/**
 * calculate phipsi codes for a peptide chain represented by the bsinchain vector
 */
std::vector<PhiPsiCodes> makephipsicodes(const std::vector<NSPgeometry::XYZ> &crds,
		const std::vector<BSInChain> & bsinchain);

/**
 * estimate SS types from phipsi
 */
std::vector<SSCode> estimatess(const std::vector<PhiPsiCodes> &phipsicodes);

/**
 * determines whether two backbone position is in the same SS segment (same helix or
 * the same beta strand)
 */
bool samesssegment(const std::vector<SSCode> &codes, int p1,int p2);
}


#endif /* IBLOCK_PHIPSICODES_H_ */

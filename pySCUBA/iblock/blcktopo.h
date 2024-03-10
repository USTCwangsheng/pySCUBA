/*
 * blcktopo.h
 *
 *  Created on: 2019年12月1日
 *      Author: hyiu
 */
#ifndef IBLOCK_BLCKTOPO_H_
#define IBLOCK_BLCKTOPO_H_
#include "geometry/xyz.h"
#include "proteinrep/aaconformer.h"
#include <vector>
#include <array>
#include <string>
#include <map>
#include <cassert>
#include <set>
#include <iostream>
#include <sstream>
namespace NSPintrct{

/**
 * Atoms are indexed using three numbers in the topology of an IntrctBlck.
 *
 * In this 3-D atom index, the first number is the chain  number,
 * the second the residue-in-chain number, and the
 * third the atom-in-residue number.
 * In BlckTopo, the chain number and residue-in-chain number are of relative meaning.
 * 0 means the current chain(residue), -1 means the previous chain(residue),
 * 1 means the next chain(residue), and so on.
 * For the atom number of an atom not in the current residue,
 * 0 correspond to N, 1 correspond to CA, 2 correspond to the first side chain atom (CB),
 * while -1 correspond to O, -2 correspond to C.
 * Note: if the index is smaller than -2, then its absolute value
 * is used to indicate the relative atom index in the other block, e.g.,
 * -3 correspond to CD1 (relative index +3),and so on.
 */
typedef std::array<int,3> AtomIdx3D;

/**
 * a simple function to help printing
 */
template <int NBODY> std::string tostring(const std::array<AtomIdx3D,NBODY> &idx){
	std::string res;
	std::ostringstream os(res);
	for(auto &i:idx) for(auto ei:i) os<<" "<<ei;
	return os.str();
}

/**
 * Bonded interactions (bonds, angles, and improper dihedrals).
 *
 * Atoms involved are indicated using the 3-D relative atom indices.
 */
template <int NBODY>
struct BondedIntrct{
	double v0; //*<lowest energy geometry value
	double k;  //*<force constant
	std::array<AtomIdx3D,NBODY> aidx; //*<participating atoms

	/**
	 * specify all data members at construction
	 */
	BondedIntrct(std::array<AtomIdx3D,NBODY> idx,double v,double kin): v0(v),k(kin),aidx(idx){;}

	/**
	 * energy calculation template method, only the specialized versions will be implemented.
	 *
	 * @param x input atomic coordinates
	 * @param dedx output energy derivatives
	 */
	double energy(const std::array<NSPgeometry::XYZ,NBODY> &x,
			std::array<NSPgeometry::XYZ,NBODY> *dedx) const;

	/**
	 * helps printing
	 */
	std::string tostring(){
		std::string res=NSPintrct::tostring<NBODY>(aidx);
		std::ostringstream oss(res,std::ostringstream::ate);
		oss <<" "<<v0<<" "<<k;
		return oss.str();
	}
};

/**
 * two atoms form a bond
 */
typedef BondedIntrct<2> Bond;

/**
 * three atoms form an angle
 */
typedef BondedIntrct<3> Angle;

/**
 * four atoms form an improper dihedral/torsion angle
 */
typedef BondedIntrct<4> ImpDih;
typedef std::array<AtomIdx3D,4> Torsion;


/**
 * a collection of various atomic properties for energy calculations and so on.
 *
 */
struct  AtomIP{

	/**
	 * calculate the steric interaction between two atoms, uses only the repulsive
	 * region of the lennard-Jones potential, with the energy minimum shifted to 0.
	 *
	 */
	static double pairenergy(const AtomIP &a1,const NSPgeometry::XYZ &x1,
			const AtomIP&a2, const NSPgeometry::XYZ &x2, std::array<NSPgeometry::XYZ,2>*dedx);
	double sigma; //*<L-J parameter, used when calculating mainchain-mainchain steric interactions
	double eps; //*<L-J parameter
	double sigmahb; //*<L-J parameter between hydrogen bond donor-acceptor atoms.
	int sterictype; //*< atom type used when calculating the sidechain packing energy term
	bool hbdonor{false};
	bool hbacceptor{false};
	bool ismainchain{false};
	bool issidechain{false};
	bool isprotein{false};
    AtomIP(int stype,double sm,double e,double shb,bool hbd,bool hba, bool mc, bool sc, bool pro):
    	sterictype(stype),sigma(sm),eps(e),sigmahb(shb),hbdonor(hbd),hbacceptor(hba),
		ismainchain(mc),issidechain(sc),isprotein(pro){;}
    std::string tostring(){
    	std::string res;
    	std::ostringstream oss(res);
    	oss<<"sterictype:"<<sterictype<<",sigma:"<<sigma<<",sigmahb:"<<sigmahb
    			<<",eps:"<<eps<<",hbdonor:"<<hbdonor<<",hbacceptor:"<<hbacceptor
				<<",ismainchain:"<<ismainchain<<",issidechain:"<<issidechain<<",isprotein:"<<isprotein;
    	return oss.str();
    }
};

/**
 * Describes atoms and covalent interactions in an IntrctBlck (e.g. an amino acid residue)
 *
 * Each BlckTopo object corresponds to one residue type. Different residues of the same type
 * can share on BlckTopo object.
 */
struct BlckTopo{
	enum BLCKPOSI{NTEMR,CTEMR,MIDDLE};//*<NTERM corresponds to n-terminal, and so on...
	static const int CD_PRO{4}; //*< the atom number index of proline CD atom
	std::string blckname; //*< usually the 3-letter residue type name.
	/**
	 * get a reference to the static BlckTopo object corresponding to the given residue type
	 *
	 * On first call, it reads a file defining the topologies and and construct the objects
	 * On subsequent calls, it returns pre-constructed objects.
	 */
	static const BlckTopo & getblcktopo_std(const std::string &residuename);

	/**
	 * Define a BlckTopo object based on a AAconformer (perhaps a small molecule ligand)
	 *
	 * Bonded terms (atoms and equilibrium geometries) are derived based on the atomic coordinates
	 * so that the molecule geometry can be closely maintained in SD simulations
	 * Atomic non-bonded parameters are given reasonable default values.
	 * Hydrogen bond donor/acceptor atom types may be specified in a separate input file.
	 */
    static const BlckTopo &getblcktopo_nonstd(const std::string &ligandname,
    		const NSPproteinrep::AAConformer *ligandconformer=nullptr);
	int natoms; //*<total number of atoms, both mainchain and sidechain
	std::vector<std::string> atomnames; //*<lists of atomnames, in the order N,CA,CB,...,C,O for amino acids
	std::map<std::string,int> atomseqs; //*<stores atom numbers indexed by atomnames
	std::vector<Bond> bonds;
	std::vector<Angle>angles;
   /**
    *  improper dihedral angles made mutable to accomodate both
    *  cis and trans peptide configuration for the same residue type.
    *  MAY NEED SPECIAL ATTENTION in future parallelization
    */
	mutable std::vector<ImpDih>impdihs;

	/**
	 *indicies of the peptide improper impdihedral angles in the
	 *impdihs vector. The lowest energy geometry of these dihedral angles
	 *impdihs varies between the cis and trans peptide bond configurations
	 */
	int pepimpdih1,pepimpdih2;

	/**
	 * torsions are used for determining 1-4 atom pairs,
	 * which are excluded from L-J steric calculations
	 */
	std::vector<Torsion> torsions;
	std::vector<AtomIP> atomips; //*< ordered as atomnames
	std::vector<std::set<AtomIdx3D>> excld; //*<list higher-index atoms to be excluded from L-J calculations.


	/**
	 *atom pairs after excluding 1-2,1-3, 1-4 pairs within the same residue
	 */
	std::vector<std::set<int>> nblist_internal_M; //*<mainchain atoms only
	std::vector<std::set<int>> nblist_internal_S; //*< involves side chain atoms
	int atomindex(const std::string  & aname) const {
		auto it=std::find(atomnames.cbegin(),atomnames.cend(),aname);
		if(it!=atomnames.cend()) return it-atomnames.cbegin();
		return -1;
	}
	std::string tostring(){
		std::string res;
		std::ostringstream oss(res);
		oss<<blckname<<" "<<natoms<<std::endl;
	    for(int i=0;i<natoms;++i){
	    	oss<<i<<" "<<atomnames[i]<<" "<<atomips[i].tostring();
	    	for(auto &x:excld[i]) for(auto ei:x) oss<<" "<<ei;
	    	oss<<std::endl;
	    }
	    for(auto &b:bonds) oss<<std::endl<<b.tostring();
	    for(auto &b:angles) oss<<std::endl<<b.tostring();
	    for(auto &b:impdihs) oss<<std::endl<<b.tostring();
	    for(auto &b:torsions) oss<<std::endl<<NSPintrct::tostring<4>(b);
	    return oss.str();
	  }
};

inline  bool less(const AtomIdx3D & i3,const AtomIdx3D & j3){
	if(i3[0]<j3[0]) return true;
	else if(i3[0]>j3[0]) return false;
	else if(i3[1]<j3[1]) return true;
	else if(i3[1]>j3[1]) return false;
	else if(i3[2]<j3[2]) return true;
	return false;
}

/**
 * add atoms to excluded list based on a bonded interaction
 */
template<int NBODY> void setexcld(std::vector<std::set<AtomIdx3D>> &excld, const std::array<AtomIdx3D,NBODY> & bd){
	for(int i=0;i<NBODY-1;++i){
		AtomIdx3D i3d=bd[i];
		for(int j=i+1;j<NBODY;++j){
			AtomIdx3D j3d=bd[j];
			if((i3d[1]!=0) &&(j3d[1]!=0)) continue; //*<at least one atom should belong to current block
			if(less(i3d,j3d)) excld[i3d[2]].insert(j3d);
			else excld[j3d[2]].insert(i3d);
		}
	}
}
}
#endif /* IBLOCK_BLCKTOPO_H_ */

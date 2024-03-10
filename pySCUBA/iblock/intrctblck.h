/*
 * iblock.h
 *
 *  Created on: 2019年12月1日
 *      Author: hyiu
 */

#ifndef IBLOCK_INTRCTBLCK_H_
#define IBLOCK_INTRCTBLCK_H_
#include "dstl/vec2d.h"
#include "geometry/xyz.h"
#include "dataio/parameters.h"
#include "iblock/blcktopo.h"
#include "designseq/ResName.h"
#include "iblock/intrctparam.h"
namespace NSPproteinrep{
class PdbRecord;
}
namespace NSPintrct{

struct IntrctBlck;

/**
 * A molecular system object is a collection of interacting block objects.
 *
 * This is the major class holding data about structure and interactions of the molecule system.
 * It is a wrapped container (vector of vectors of IntrctBlck  objects).
 * In the case of proteins, each IntrctBlck object is an amino acid residue,
 * each vector of intrctblck objects is a peptide chain,and the entire molecule comprises
 * one or more chains.
 */
typedef NSPdstl::Vec2D<IntrctBlck> MolSystm;

/**
 * Extra bond angles and improper dihedral angles at a pre-PRO peptide bond
 */
struct PrePro{
   std::vector<Angle>angles; //*<extra bond angles in e.g. preproline
   std::vector<ImpDih>impdihs; //*<extra impdihs in e.g. preproline
   PrePro(const BlckTopo &);
   static void extraforces(const IntrctBlck &blck, const IntrctPara &para,
		   std::vector<NSPgeometry::XYZ> *dedx_sys);
};

class IntrctMol;

/**
 * A molecular system is considered as a collection of interaction blocks.
 *
 * A block usually correspond to an amino acid residue within a peptide bond
 * It may also hold a small molecule in a complex
 * Bonded interactions across blocks are included in blcktopo using
 * 3-D relative atom indices. The absolute atom indices are deduced from
 * the relative atom indices and block-specific offsets
 */
struct IntrctBlck{
public:
	enum ENES{BOND,ANGLE,IMPDIH,STERIC,PHIPSI,LOCALSTR,LOCALHB,SITEPAIR,ROTAMER,
		SCPACKING,ENESIZE}; //<* locations of different energy components in the energy array
	enum ACTIVE{ALL,SIDECHAIN,INACTIVE}; //<*all atom can move, only side chain atoms can move or all fixed in sampling
	enum REGION{PROTEIN_POCKET=-1,PROTEIN_CORE,PROTEIN_VAR,LIGAND}; //<* PROTEIN_VAR will not see PROTEIN_POCKET in pocket design simulations
	const BlckTopo *topo{nullptr};  //*<topo data may be shared between different residues of the same  type
	int aoffset{0}; //*<for deducing absolute atom indices
    int activemod{ALL}; //*<indicates how coordinates changes in confomation sampling
    int regioncode{PROTEIN_CORE};//*< blocks of regioncode <0 will not see blocks of regioncode >0 for non-bonded forces
    std::set<int> lexcludedatoms;//*<protein atoms that are excluded from steric interactions with any ligand atom;
    double weight_mc{1.0}; //*<weight of the block for mainchain local HB and sitepair  energy calculation
    mutable std::vector<NSPgeometry::XYZ> crds; //*<atomic coordinates, atom orders are defined in topo
    bool isaminoacid{true}; //*<flag
    bool precis{false}; //*< flags if the residue forms a cis peptide bond with next residue;
    bool prepro{false}; //*< flags if the next residue is a proline;
    bool softsc{false}; //*< flags if the side chain packing is estimated using a soften function
	mutable std::array<double,ENESIZE> energies; //*<stores energy components computed by the various force calculating methods
	typedef std::vector<std::set<int>> AtomNbList;   //*< stores for each atom of this block neighboring atoms in another block.
	typedef std::map<NSPdstl::Idx2D,AtomNbList> BlckAtomNbList; //*<Neighboring atoms in neighboring blocks.
    mutable BlckAtomNbList nblist_M,nblist_S,nblist_L; //*< _M: between mainchain atoms only, _S: at least one side chain atom
    																											//*<_L: involving non-protein(ligand) atoms
    IntrctBlck(){;}

    void insert_nblist_M(const IntrctBlck &blck2, const AtomNbList &list) const {
		nblist_M[blck2.getidx2d()]=list;
	}
	void insert_nblist_S(const IntrctBlck &blck2, const AtomNbList &list) const {
		nblist_S[blck2.getidx2d()]=list;
	}
	void insert_nblist_L(const IntrctBlck &blck2, const AtomNbList &list) const {
		nblist_L[blck2.getidx2d()]=list;
	}

/**
 * obtain block index offset of the chain containing the block
 *
 * The absolute block index is blkoffset() + idx2d_.first
 */
	int blkoffset() const{
		if(molsystm_) return molsystm_->offset(idx2d_.first);
		else return 0;
	}

	/**
	 * set pointer to the containing molecule, the block's own chain index, and its
	 * relative index in the chain
	 */
	void setimol(IntrctMol * mol,int id1,int id2);

	/**
	 * set atomic coordinates by copying from a source
	 */
	void copycrdsfrom(std::vector<NSPgeometry::XYZ>::const_iterator it) const {
		auto it_end=it+crds.size();
		std::copy(it,it_end,crds.begin());
	}
	/**
	 * copy atomic coordinates to a destination
	 */
	void copycrdsto(std::vector<NSPgeometry::XYZ>::iterator it) const {
		std::copy(crds.begin(),crds.end(),it);
	}

	/**
	 * copy atomnames
	 */
	void copyatomnamesto(std::vector<std::string>::iterator it) const{
		std::copy(topo->atomnames.begin(),topo->atomnames.end(),it);
	}

	/**
	 * get pointer to the containing molecule object
	 */
	IntrctMol *getmol() const {return imol_;}

	/**
	 * get pointer to the container molecular object
	 */
	MolSystm *getsystm() const {return molsystm_;}

	/**
	 * idx2d.first is the chain index,  idx2d.second is the relative index of the block in chain
	 */
	NSPdstl::Idx2D getidx2d() const{return idx2d_;}

	/**
	 * get the chain index
	 */
	int chainseq() const{return idx2d_.first;}

	/**
	 * get relative index in chain
	 */
	int resseq() const {return idx2d_.second;}

	/**
	 * get the name of atom i
	 */
	std::string atomname(int i) const {return topo->atomnames.at(i);}

	/**
	 * get the coordinates of atom i
	 */
	NSPgeometry::XYZ getcrd(int i) const {return crds.at(i);}

	/**
	 * get atom coordinates by atom name
	 */
	NSPgeometry::XYZ getcrd(const std::string &atomname) const {
		return crds[topo->atomseqs.at(atomname)];
	}

	/**
	 * get a coordinate point to represent the block.
	 *
	 * If the block is a amino acid residue, it's the coordinate of the CA atom.
	 * Otherwise it's the geometric center of all atoms
	 */
	NSPgeometry::XYZ getrepcrd()const ;
	IntrctBlck *nextblck() const{
		if(molsystm_) {
			if(idx2d_.second>=molsystm_->sized2(idx2d_.second)-1) return nullptr;
			return & (*molsystm_)(idx2d_.first,idx2d_.second+1);
		}
		return nullptr;
	}
	IntrctBlck *prevblck() const{
		if(molsystm_ && idx2d_.second>0) return & (*molsystm_)(idx2d_.first,idx2d_.second-1);
		return nullptr;
	}

	/**
	 * number of atoms in this block, including both main chain and side chain atoms
	 */
	int natoms()const {if(topo) return topo->natoms; else return 0;}

	/**
	 * residue name
	 */
	std::string resname() const {if(topo) return topo->blckname;return "";}

	/**
	 * absolute index of backbone atoms (N,CA,C,O) in this block
	 * determined using aoffset and relative atom index
	 */
	std::vector<int> getbackboneatoms() const{
		if(!NSPdesignseq::ResName::resname().isStandardAminoAcid(resname())) return std::vector<int>();
		return std::vector<int>{aoffset,aoffset+1,aoffset+natoms()-2,aoffset+natoms()-1};
	}
	/**
	 * Make a backbone site object from the block.
	 *
	 *Only the residue name and coordinates will be set. Other data
	 *such as phipsi angles are NOT set here.
	 */
	NSPproteinrep::BackBoneSite getbackbonesite() const;
	/**
	 * determine atomic neighbor lists.
	 */
	void mksteric_nblist(const IntrctPara &param) const;

	/**
	 * change coordinates of side chain atoms
	 */
	void changesccrds(const std::vector<NSPgeometry::XYZ> &sccrds);

	/**
	 * printing for debug
	 */
	std::string print_nblist() const {
		std::ostringstream oss;
		for(auto &l:nblist_M){
			oss<<"block"<<l.first.first<<":"<<l.first.second<<std::endl;
			for(int i=0;i<l.second.size();++i){
				for(auto j:l.second[i]) oss<<" "<<j;
				oss<<std::endl;
			}
		}
		return oss.str();
	}

	/**
	 * determine the absolute atom index by given a 3-D relative index
	 */
	int getatomindex_abs(const AtomIdx3D & j3d) const {
		NSPdstl::Idx2D bidx2d{idx2d_.first+j3d[0],idx2d_.second+j3d[1]}; //change  relative 2d index to absolute 2d index
		int joffset=(*molsystm_)(bidx2d.first,bidx2d.second).aoffset;
		int res;
		if(j3d[2]>=0) res=joffset+j3d[2];
		else if(j3d[2]>=-2)res=joffset+(*molsystm_)[bidx2d].natoms()+j3d[2];
		else res =joffset-j3d[2];
		return res;
	}
	/**
	 * get a 3-D absolute index from a relative 3D index
	 */
	AtomIdx3D toabs_index(const AtomIdx3D &j3d) const{
		int i1=idx2d_.first+j3d[0],i2=idx2d_.second+j3d[1],i3=j3d[2];
		if(i3<-2) i3=-i3;
		else if(i3<0) i3=(*molsystm_)(i1,i2).natoms()+i3;
		return AtomIdx3D{i1,i2,i3};
	}
	/**
	 * get a 3-D relative index from an absolute one
	 */
	AtomIdx3D torelative_index(const AtomIdx3D &j3d) const{
		int natomsb=(*molsystm_)(j3d[0],j3d[1]).natoms();
		int i1=j3d[0]-idx2d_.first,i2=j3d[1]-idx2d_.second,i3=j3d[2];
		if(i1==0 && i2==1) {
			if(i3>1) {
				if(i3<natomsb-2) i3=-i3;
				else i3=i3-natomsb;
			}
		}
		return AtomIdx3D{i1,i2,i3};
	}

	/**
	 * helper function for  printing
	 */
	std::pair<std::string,std::string> getnames(const AtomIdx3D &j3d) const{
		NSPdstl::Idx2D bidx2d{idx2d_.first+j3d[0],idx2d_.second+j3d[1]};
		return std::make_pair((*molsystm_)[bidx2d].resname(),(*molsystm_)[bidx2d].atomname(j3d[2]));
	}
	/**
	 * printing helper
	 */
	std::string atomstring(int i) const{
		std::ostringstream oss;
		oss<<"chain"<<idx2d_.first<<resname()<<idx2d_.second<<atomname(i);
		return oss.str();
	}

	/**
	 * generate PdbRecords for atoms in the block, for writing a PDB file and so on
	 */
	std::vector<NSPproteinrep::PdbRecord> topdbrecords()const;

	/**
	 * calculate a group of bonded interactions of the same class
	 */
	template<int NBODY> double bondedforce(const IntrctPara &param,const std::vector<BondedIntrct<NBODY>> & terms,
			std::vector<NSPgeometry::XYZ> *dedx_sys) const{
		double esum=0.0;
		double w=param.weight_coval;
		std::ofstream *ofs{nullptr};
		if(param.enedetails)ofs=&(Results::ofstreams(param.jobname,"CovalentEne"));
		for(auto & term:terms){
			std::array<NSPgeometry::XYZ,NBODY> r,dedx;
			bool exist=true;
			for(int i=0;i<NBODY;++i) {
				AtomIdx3D i3d=term.aidx.at(i);
//				if(i3d[0]==0 && i3d[1]==0){
//					r[i]=crds[term.aidx.at(i)[2]];
//				} else {
					NSPdstl::Idx2D bidx2d{idx2d_.first+i3d[0],idx2d_.second+i3d[1]};
					if(bidx2d.second>=molsystm_->sized2(bidx2d.first)) {
						exist=false;
						break;
					}
					int i3=i3d[2];
					if(i3<-2) i3=-i3;
					else if(i3<0)i3=(*molsystm_)[bidx2d].natoms()+i3;
					r[i]=(*molsystm_)[bidx2d].getcrd(i3);
//				}
			}
			if(!exist) break;
			double e=w*term.energy(r,&dedx);

			esum+=e;
			for(int i=0;i<NBODY;++i){
		    	int iidx= getatomindex_abs(term.aidx.at(i));
		        (*dedx_sys) [iidx] = (*dedx_sys)[iidx]+w*dedx[i];
		        auto names=getnames(term.aidx.at(i));
		        if(ofs)  (*ofs) <<" "<<term.aidx.at(i)[0]<<"-"<<term.aidx.at(i)[1]<<"-"<<names.first<<"-"<<names.second;
		    }
			if(ofs) (*ofs) <<" e = "<<e<<std::endl;
		}
		return esum;
	}
	/**
	 * calculate all bonded forces
	 */
	double covalentforces(const IntrctPara &param,std::vector<NSPgeometry::XYZ> *dedx_sys) const {
		energies[BOND]=bondedforce(param,topo->bonds,dedx_sys);
		energies[ANGLE]=bondedforce(param,topo->angles,dedx_sys);
		if(precis){
			topo->impdihs[topo->pepimpdih1].v0=0.0;
			topo->impdihs[topo->pepimpdih2].v0=3.14159265;
		}
		energies[IMPDIH]=bondedforce(param,topo->impdihs,dedx_sys);
    	if(precis){
				topo->impdihs[topo->pepimpdih1].v0=3.14159265;
				topo->impdihs[topo->pepimpdih2].v0=0.0;
			}
		if(prepro){
			PrePro::extraforces(*this,param,dedx_sys);
		}
		return energies[BOND]+energies[ANGLE]+energies[IMPDIH];
	}

	/**
	 * calculate mainchain-mainchain steric forces of the block
	 * with another block in the molecular system
	 */
	double blckstericforces(const IntrctPara &param, const IntrctBlck &blk2,const AtomNbList &nblist,
			std::vector<NSPgeometry::XYZ> *dedx_sys) const;

	/**
	 * calculate sidechain packing forces (sidechain-mainchain and sidechain-sidechain)
	 * with another block in the molecular system
	 */
	double blkscpackingforces(const IntrctPara &param,
			const IntrctBlck &blk2, const AtomNbList &nblist,
			std::vector<NSPgeometry::XYZ> *dedx_sys) const;

	/**
	 * calculate mainchain-mainchain steric forces with all subsequent blocks in the
	 * molecular system.
	 *
	 * It calls blckstericforces
	 */
    double stericforces(const IntrctPara &param, std::vector<NSPgeometry::XYZ> *dedx_sys) const;

    /**
     * calculate sidechain-mainchain and sidechain-sidechain packing forces with
     * all subsequent blocks in the molecular system.
     *
     * It calls blkscpackingforces
     */
    double scpackingforces(const IntrctPara &param, std::vector<NSPgeometry::XYZ> *dedx_sys) const;

    /**
     * determines atomic neighbors from another block
     * The first element of the returned vector contains mainchain-mainchain atomic neighbors
     * The second element contains sidechain-sidechain and sidechain-mainchain (or mainchain-sidechain)
     * neighbors
     * The third element contains protein-nonprotein neighbors
     * nnb delivers the number of neighbors of each type.
     */
	std::vector<AtomNbList> getneighbors(const IntrctBlck &blck2,
	   const IntrctPara &param,std::array<int,3> &nnb) const;
private:
    IntrctMol *imol_{nullptr}; //*< points to the containing the IntrctMol object
	MolSystm * molsystm_{nullptr}; //*< points to the containing MolSystm container. redundant wrt. imol_
	NSPdstl::Idx2D idx2d_{-1,-1}; //*< 2D index of self in molsystm, firt=chainid, second=relative residue id
};
}




#endif /* IBLOCK_INTRCTBLCK_H_ */

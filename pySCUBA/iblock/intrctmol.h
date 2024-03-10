/*
 * intrctmol.h
 *
 *  Created on: 2019年12月3日
 *      Author: hyiu
 */

#ifndef IBLOCK_INTRCTMOL_H_
#define IBLOCK_INTRCTMOL_H_
#include "iblock/intrctblck.h"
#include "iblock/aasequtil.h"
//#include "fullsite/fullsite.h"
#include "proteinrep/aaconformer.h"
#include "iblock/phipsicodes.h"
#include "iblock/vsctypes.h"
#include "iblock/blkselector.h"
#include "include/debuglevel.h"
#define DLEVEL_CHECKINDEX 7
#define DLEVEL_CHECKPEPTIDE 7
#define DLEVEL_CHECKCRDS 9
namespace  NSPintrct{
/**
 * The main class describes structure and interaction of a molecular system.
 *
 * Many member objects are organized as vector of vectors of residues.
 * The outer layer elements refer to peptide chains.
 * The inner layer elements refer to residues in a chain.
 * The program testintrctmol illustrates how to use this class
 */
class IntrctMol{
public:
	/**
	 * peptide backbone atoms in the molecule
	 *
	 * We do not assume all chains correspond to peptides.
	 */
	std::vector<std::vector<BSInChain>> bsinchains;

	/**
	 * sidechain atoms in the molecule
	 */
	std::vector<std::vector<SCInChain>> scinchains;

	/**
	 * phipsi and encoded values calculated on a configuration.
	 */
	mutable std::vector<std::vector<PhiPsiCodes>> phipsicodes;

	/**
	 * secondary structures and encoded values calculated on a configuration
	 */
	mutable std::vector<std::vector<SSCode>> sscodes;

	/**
	 * sidechain torsions and encoded values calculated on a configuration
	 */
	mutable std::vector<std::vector<ConformerCode>> conformercodes;

	/**
	 * internal control parameters for force calculation
	 *
	 * In some method, externally provided control parameters are used instead of
	 * these internal parameters. In other words, the externally provided
	 * parameters are of higher priority than the internal ones.
	 */
	mutable std::shared_ptr<IntrctPara> my_intrctpara{nullptr};

	IntrctMol(){;}

	/**
	 * copy constructor has been changed, to update the container-pointing pointers in
	 * the contained intrctblck objects
	 */
	IntrctMol(const IntrctMol &mola);

	/**
	 * set up using amino acid sequences represented by 1-letter code strings.
	 * Each string represents a peptide chain.
	 * Coordinates will be left undefined.
	 */
	IntrctMol(const std::vector<std::string> chains1l){
		setup(chains1l);
	}

	/**
	 * setup using amino acid sequences represented by a vector of vectors of 3-letter
	 * residue type string.
	 * Each vector of strings represent a peptide chain
	 * Coordinates will be left undefined.
	 *
	 */
	IntrctMol(const std::vector<std::vector<std::string>> chains3l){
		setup(chains3l);
	}
/*
	/**
	 * set up using a vector of vector of FullSite objects.
	 *
	 * Each fullsite object correspond to a residue (e.g. read in from a PDB file)
	 * Residue names and atomic coordinates in the full site objects will be used.
	 *
	 * Note:
	 * if a fullsite object does not contain all needed sidechain
	 * atom coordinates, the corresponding residue will be changed into a glycine.
	 * This behavior may need to be changed later.

	IntrctMol(const std::vector<std::vector<NSPproteinrep::FullSite>> &chains){
		setup(chains);
	}
*/
	/**
		 * set up using a vector of vector of AAConformer objects.
		 *
		 * Each conformer object correspond to a residue (e.g. read in from a PDB file)
		 * Residue names and atomic coordinates in the full site objects will be used.
		 *
		 * Note:
		 * if a conformer object does not contain all needed sidechain
		 * atom coordinates, the corresponding residue will be changed into a glycine.
		 * This behavior may need to be changed later.
		 */
		IntrctMol(const std::vector<std::vector<NSPproteinrep::AAConformer>> &chains,
				bool ignore_nonprotein=true){
			setup(chains,ignore_nonprotein);
		}

		IntrctMol(const NSPproteinrep::AAConformersInModel  &model){
				setup(model.conformers);
				this->mappdbkeyint_=model.mappdbkeyint;
			}

	/**
	 * the actual setup methods called by the constructors
	 */
	void setup(const std::vector<std::vector<std::string>> & chains3l);

	/**
		 * the actual setup methods called by the constructors
		 */
	void setup(const std::vector<std::vector<NSPproteinrep::AAConformer>> &chains,bool ignore_nonprotein=true);
//	void setup(const std::vector<std::vector<NSPproteinrep::FullSite>> &chains);

	/**
		 * the actual setup methods called by the constructors
		 */
	void setup(const std::vector<std::string> & chains1l){
		std::vector<std::vector<std::string>> c3l;
		for(int c=0;c<chains1l.size();++c){
		    c3l.push_back(getseq3l(chains1l.at(c)));
		}
		setup(c3l);
	}
	void addchain(const std::vector<IntrctBlck> &chain){
		mol_.push_backvec(chain);
		auto & c=mol_.backvec();
		int cidx=mol_.sized1()-1;
		for(int i=0;i<c.size();++i){
			c[i].setimol(this,cidx,i);
		}
		this->offsetsok_=false;
	}
	/**
	 * replace the residue at a position by another one of another type.
	 * The topo is changed. Coordinates not setup in this method.
	 */
	void replaceblck(NSPdstl::Idx2D idx2d,const std::string &residuename);

	/**
	 * replace the intrctblck by another one (which may contain new coordinates)
	 *
	 * Used for mutation and so on.
	 */
	void pasteblck(const IntrctBlck &blk);

	/**
	 * change coordinates of all atoms.
	 *
	 * To use the const version, the atom index offsets of the blocks must have already
	 * been correctly set.
	 */
	void changecrds(const std::vector<NSPgeometry::XYZ> &newcrds) const;

	/**
	 * change coordinates of all atoms.
	 *
	 * This non-const version can set the atom index offsets if needed.
	 */
	void changecrds(const std::vector<NSPgeometry::XYZ> &newcrds) {
		if(!offsetsok_) setblckoffsets(NSPdstl::Idx2D{0,0});
		((const IntrctMol*) this)->changecrds(newcrds);
	}

	MolSystm & getmolsystm() {return mol_;}
	const MolSystm &getmolsystm() const{return mol_;}

	/**
	 * number of chains
	 */
	int nchains() const {return mol_.sized1();}

	/**
	 * number of blocks(residues) in a chain
	 */
	int nresidues(int c) const{return mol_.sized2(c);}

	/**
	 * total number of atoms.
	 *
	 * To use this const version, the atom index offsets of all blocks must have
	 * already been set properly.
	 */
	int natoms() const{
#if DEBUG_LEVEL > DLEVEL_CHECKINDEX
			assert(checkblckindices());
			assert(checkaoffsets());
#endif
		assert(offsetsok_);
		return mol_.backvec().back().aoffset+mol_.backvec().back().natoms();}

	/**
	 * total number of atoms.
	 *
	 * This non-const version can set the atom index off blocks if needed.
	 */
	int natoms(){
		if(offsetsok_) {
#if DEBUG_LEVEL > DLEVEL_CHECKINDEX
			assert(checkblckindices());
			assert(checkaoffsets());
#endif
			return mol_.backvec().back().aoffset+mol_.backvec().back().natoms();
		}
		setblckoffsets(NSPdstl::Idx2D{0,0});
		return mol_.backvec().back().aoffset+mol_.backvec().back().natoms();
	}

	IntrctBlck & getblck(const NSPdstl::Idx2D &idx2d){
		return mol_[idx2d];
	}
	const IntrctBlck & getblck(const NSPdstl::Idx2D &idx2d) const{
		return mol_[idx2d];
	}
	IntrctBlck & operator()(const NSPdstl::Idx2D &idx2d){
		return mol_[idx2d];
	}
	const IntrctBlck & operator()(const NSPdstl::Idx2D &idx2d)const {
			return mol_[idx2d];
		}

	/**
	 * get atomname
	 * i is the absolute 1-D atom index
	 */
	std::string atomname(int i) const {
		const IntrctBlck &blk=(*this)(residueofatoms_[i]);
		return blk.atomname(i-blk.aoffset);
	}

	std::string residuename(NSPdstl::Idx2D ridx) const {
		return mol_[ridx].resname();
	}

	/**
	 * returns a vector contains the 2-D residue index of all atoms
	 *
	 * To use this constant version, the member data should have already been properly set.
	 */
	const std::vector<NSPdstl::Idx2D> &residueofatoms() const{
#if DEBUG_LEVEL > DLEVEL_CHECKINDEX
			assert(checkblckindices());
			assert(checkaoffsets());
#endif
		assert(offsetsok_);
		return residueofatoms_;
	}

	/**
	 * Returns the 2-D residue indices of all atoms
	 *
	 * This non-const version can setup the atom index offsets if needed
	 *
	 */
	const std::vector<NSPdstl::Idx2D> &residueofatoms(){
		if(offsetsok_) {
#if DEBUG_LEVEL > DLEVEL_CHECKINDEX
			assert(checkblckindices());
			assert(checkaoffsets());
#endif
			return residueofatoms_;
		}
		setblckoffsets(NSPdstl::Idx2D{0,0});
		return residueofatoms_;
	}
    /**
     * get a vector of atom indexes contained a list of given blocks
     */
	std::vector<int> atomindices(const std::map<int,std::set<int>> & chainblcks) const;

 /**
	 * get the coordinates of all atoms as a vector of XYZ objects.
	 *
	 * This const version requires that the atom offsets of blocks have already been set properly.
	 */
   const std::vector<NSPgeometry::XYZ> & recollectcrds_all() const;

  /**
   * translate a chain of chainid so that the geometry center of all atoms is at  newcenter
   */
   void translatechaincenter(const NSPgeometry::XYZ &newcenter, int chainid=0);

   /**
    * get the coordinates of all atoms as a vector of XYZ objects.
    *
    * This non-const version can setup the atom offsets if needed.
    */
	const std::vector<NSPgeometry::XYZ> & recollectcrds_all() {
		if(!offsetsok_) setblckoffsets(NSPdstl::Idx2D{0,0});
		return ((const IntrctMol*) this)->recollectcrds_all();
	}

	/**
	 * get the coordinates of all atoms as vector of doubles (size=3*natoms)
	 *
	 * const version requires atom offsets to have been properly set.
	 */
	const std::vector<double> & getcrds_all_d() const{
		recollectcrds_all();
		return crds_all_d_;
	}

	/**
	 * get the coordinates of all atoms as a vector of doubles (size=3*natoms)
	 *
	 * non-const version can set the atom offsets if required.
	 */
	const std::vector<double> & getcrds_all_d(){
		recollectcrds_all();
		return crds_all_d_;
	}

	/**
	 * returns a vector of bools flagging whether each an atom moves in a sampling process.
	 *
	 * The bool values are determined based on the activemod of each block.
	 * The difference between the const and the non-const version is the same as the other members.
	 */
	const std::vector<bool> & isfixedvec() const{
#if DEBUG_LEVEL > DLEVEL_CHECKINDEX
			assert(checkblckindices());
			assert(checkaoffsets());
			assert(checkatomsfixed());
#endif
		assert(offsetsok_);
		assert(isatomfixed_.size()==natoms());
		return isatomfixed_;
	}

	/**
		 * returns a vector of bools flagging whether each an atom moves in a sampling process.
		 *
		 * The bool values are determined based on the activemod of each block.
		 * The difference between the const and the non-const version is the same as the other members.
		 */
	const std::vector<bool> & isfixedvec(){
		if(!offsetsok_) setblckoffsets({0,0});
		if (isatomfixed_.size()==natoms()) {
#if DEBUG_LEVEL > DLEVEL_CHECKINDEX
			assert(checkatomsfixed());
#endif
			return isatomfixed_;
		}
		return setatomfixed();
	}
	/**
	 * set the activemod of all blcks to mode
	 */
	void changeactivemodeveryblck(int mode){
		for(auto &c:mol_.D_){
			for(auto &b:c){
				b.activemod=mode;
			}
		}
		offsetsok_=false;
	}
	void fixeveryblck(){
		changeactivemodeveryblck(IntrctBlck::INACTIVE);
	}
	void activateeveryblck(){
		changeactivemodeveryblck(IntrctBlck::ALL);
	}
	/**
	 * Secify activemod of IntrctBlck objects of this mol.
	 *
	 * First, set all blcks as fixed. Then
	 * set the activemod of selected blcks in mcsel to allactive
	 * and those in scsel but not in mcsel to sidechain active.
	 */
	void specifyactiveblks(const BlkSelector::SelectedBlks &mcsel,
			const BlkSelector::SelectedBlks &scsel){
		this->fixeveryblck();
		for(auto &aset:scsel){
			for(int r:aset.second){
				(*this)({aset.first,r}).activemod=IntrctBlck::SIDECHAIN;
			}
		}
		for(auto &aset:mcsel){
			for(int r:aset.second){
				(*this)({aset.first,r}).activemod=IntrctBlck::ALL;
			}
		}
	}
	/**
	 * Specify activemod of IntrctBlck objects of this mol.
	 *
	 * First, set all blcks as all active. Then
	 * set the activemod of selected blcks in allfixedsel to inactive
	 * and those in mcfixedsel but not in allfixedsel to sidechain active.
	 */
	void specifyfixedblks(const BlkSelector::SelectedBlks &allfixedsel,
			const BlkSelector::SelectedBlks &mcfixedsel){
		this->activateeveryblck();
		for(auto &aset:mcfixedsel){
			for(int r:aset.second){
				(*this)({aset.first,r}).activemod=IntrctBlck::SIDECHAIN;
			}
		}
		for(auto &aset:allfixedsel){
			for(int r:aset.second){
				(*this)({aset.first,r}).activemod=IntrctBlck::INACTIVE;
			}
		}
	}
	void calcphipsicodes(bool updateonly=false) const;
	void calcsscodes(bool updateonly=false) const;
	std::vector<std::string> sscodestrings() const;
	void calcconformercodes(bool updateonly=false) const;
	void clear_energies();
	double sum_energies(std::array<double,IntrctBlck::ENESIZE> &energies) const;

/**
 * a helper determine if there is any main chain atom within a window around a residue position
 * can move during a simulation
 *
 */
	bool mainchainwindowactive(const NSPdstl::Idx2D &posi, int wh) const{
		int first=posi.second-wh;
		if(first<0) first=0;
		int last=posi.second+wh+1;
		if(last>=mol_.D_[posi.first].size()) last=mol_.D_[posi.first].size();
		for(int i=first;i<last;++i)
			if(mol_(posi.first,i).activemod==IntrctBlck::ALL) return true;
		return false;
	}

	/**
			 * Main entrance calculating total energy and total forces.
			 *
			 * It calls methods to calculate different energy terms, sum up the forces.
			 * The energies are stored in individual blocks, which can be summed up
			 * by a subsequent calling of the sum_energies method
			 *
			 * The supplied, not the internal interaction parameters will be used.
			 */
		void forces_all(const IntrctPara &param,
				std::vector<NSPgeometry::XYZ> &dedx) const;

	/**
	 * this overload wraps the force calculation by supplying the internal
	 * internal interaction parameter.
	 */
	void forces_all(std::vector<NSPgeometry::XYZ> &dedx) const{
		forces_all(*my_intrctpara,dedx);
	}

	/**
		 * this overload wraps the force calculation by changing the atomic coordinates
		 * before force calculations and use the internal interaction parameter.
		 */
	void forces_all(const std::vector<NSPgeometry::XYZ> &crds,
			std::vector<NSPgeometry::XYZ> &dedx) const{
		forces_all(*my_intrctpara,crds,dedx);
	}

	/**
		 * this overload wraps the force calculation by changing the atomic coordinates
		 * before the force calculation
		 */
	void forces_all(const IntrctPara &param,
			const std::vector<NSPgeometry::XYZ> &crds,
			std::vector<NSPgeometry::XYZ> &dedx) const{
		this->changecrds(crds);
		forces_all(param,dedx);
	}

	void forces_cov(const IntrctPara &param,std::vector<NSPgeometry::XYZ> &dedx) const;
	void mksteric_nblists(const IntrctPara &param) const;
	void forces_steric(const IntrctPara &param,std::vector<NSPgeometry::XYZ> &dedx) const;
	void forces_phipsi(const IntrctPara &param,std::vector<NSPgeometry::XYZ> &dedx) const;
	void forces_localstr(const IntrctPara &param,std::vector<NSPgeometry::XYZ> &dedx) const;
	void forces_sitepair(const IntrctPara &param,std::vector<NSPgeometry::XYZ> &dedx) const;
	void forces_rotamer(const IntrctPara &param, std::vector<NSPgeometry::XYZ> &dedx) const;
	void forces_scpacking(const IntrctPara &param, std::vector<NSPgeometry::XYZ> &dedx) const;
	void forces_localhb(const IntrctPara &param, std::vector<NSPgeometry::XYZ> &dedx) const;
     void ssregions(int chainid, std::vector<std::pair<int,int>> &helixregions,
    		 std::vector<std::pair<int,int>> &strandregions);
	/**
	 * write the structure to a PDB file
	 */
	void writepdb(std::ofstream &ofs) const;
/**
 * check consistency of block indices, used for debugging purposes
 */
	bool checkblckindices() const;
/*
 * check consistency of atomoffsets of blocks, used for debugging purposes
 */
	bool checkaoffsets() const;
/**
 * check bsinchains, used for debugging
 */
	bool checkbsinchains() const;
/**
 * check atomsfixed, used for debugging
 */
	bool checkatomsfixed() const;
/**
 * check scinchains, used fort debugging
 */
	bool checkscinchains() const;
/**
 * check consistency between global coordinate and block coordinate, usde for debugging
 */
	bool checkcrds() const;
	std::shared_ptr<NSPproteinrep::MapPdbKeyInt> & mappdbkeyint() {return mappdbkeyint_;}
	const std::shared_ptr<NSPproteinrep::MapPdbKeyInt> & mappdbkeyint() const  {return mappdbkeyint_;}
private:
	MolSystm mol_; //*< the wrapped Vec2D container holding all the interctblck objects.
	std::shared_ptr<NSPproteinrep::MapPdbKeyInt> mappdbkeyint_ { nullptr}; //maps for chainids and residue id in original PDB
	bool offsetsok_{false}; //*< flags whether atom index offsets of the blocks have been properly set
	std::vector<NSPdstl::Idx2D> residueofatoms_; //*<2-D residue index of all atoms
	std::vector<bool> isatomfixed_; //*<flags if an atom moves in configuration sampling based on the block activemode
	mutable std::vector<double> crds_all_d_; //*< coordinates of all atoms as vector of doubles (size=3*natoms)
	mutable std::vector<NSPgeometry::XYZ> crds_all_; //*< coordinates of all atoms as vector of XYZ objects.
	const std::vector<bool> & setatomfixed(); //*<set up the isatomfixed_ data member using the activemode of blocks
	void setblckoffsets(const NSPdstl::Idx2D & startblck); //*<setup the atom index offsets of blocks
	void localstrenergy(const IntrctPara &param,
				const std::vector<BSInChain> & bsinchain,
				const std::vector<PhiPsiCodes> & phipsicodes,
				std::vector<NSPgeometry::XYZ> *dedx_sys) const; //*<calculate the local structure energy of one chain
};
struct MolSystmPara;

/**
 * make a molecule system based on parameters and instructions in the MolSystmPara object
 * and the IntrctPara object.
 */
std::shared_ptr<IntrctMol> make_molsystm(const MolSystmPara &molpara,
		const IntrctPara &ipara);
}


#endif /* IBLOCK_INTRCTMOL_H_ */

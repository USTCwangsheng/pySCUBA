/*
 * molmodeler.h
 *
 *  Created on: 2019年12月18日
 *      Author: hyiu
 */

#ifndef IBLOCK_MOLMODELER_H_
#define IBLOCK_MOLMODELER_H_
#include "iblock/intrctmol.h"
#include "designseq/ProteinRep.h"
#include "iblock/loopreservoir.h"

namespace NSPdesignseq {
class RotamerLib;
class Rotamer;
}
namespace NSPintrct {
/**
 * For modifying the structure and sequence of peptides in a molecular system.
 *
 * The program testmolmodeler illustrates how to use this class
 */
class MolModeler {
public:
	/**
	 * 2-D index of a residue, composed of chain number and residue number.
	 */
	typedef NSPdstl::Idx2D SiteIdx;
	/**
	 * wrapper of loop reservoir
	 */
	struct LoopConfMaker{
		LoopReservoir reservoir;
		SiteIdx leftsite;
		SiteIdx rightsite;
	};

	/**
	 * set the target molecular system to modify.
	 */
	void settarget(IntrctMol &target) {
		target_ = &target;
		siteneighbors_.clear();
		phipsis_.clear();
		sitelocalframes_.clear();
		backbonesites_.clear();
	}

	/**
	 * access the backbonesites_ data
	 *
	 *This non-const overload setup the data when needed
	 */
	const std::vector<std::vector<NSPproteinrep::BackBoneSite>> &getbackbonesites();

	/**
	 * access the backbonesites_ data
	 */
	const std::vector<std::vector<NSPproteinrep::BackBoneSite>> &getbackbonesites() const{
		assert(backbonesites_.size()==target_->nchains());
		return backbonesites_;
	}

	/**
	 * set value of the checkflash_ flag
	 */
	void checkclash(bool val=true){
		checkclash_=val;
	}

	/**
	 * set the CA-CA distance cutoff for determining neighboring sites.
	 * clash check is only carried out between sites within cutoff.
	 * @param r the cutoff distance in nanometer.
	 */
	void setneigborcutoff(double r){
		neighbor_cutoff2_=r*r;
	}

	/**
	 * modify the sidechain type and conformation at a given site.
	 * A rotamer leading to the fewest atomic clashes with the rest of the molecule
	 * (if checkclash_ flag is true) and of the lowest backbone-dependent
	 * conformational energy will be used.
	 * @param site the position
	 * @resname new sidechain type
	 */
	void changeresidue(const SiteIdx &site, std::string resname){
		target_->pasteblck(buildsidechain(site,resname));
	}

	/**
	 * Copy the target and change it by replacing the loop between siteright and siteleft(
	 * not including).
	 *
	 * This method can be used when the loop length has changed
	 */
	 IntrctMol mkimol_newloop
	 	 (const SiteIdx &siteleft, const SiteIdx &siteright,
	 			 const std::vector<IntrctBlck> &newloop) const;
	 /**
	  * Add a multiple chains,
	  */
	 /**
	  * Merging two chains into one and copy the remaining chains of target_.
	  *
	  * Note: the newly merged chain IS NOT properly connected.
	  * On returning, the leftsite is delivered with the new index of
	  * the residue ahead of the C-terminal residue of the first chain,
	  *  while the rightsite the new index of the residue after the N-terminal residue
	  *  of the second chain. Later on, a linker of user-specified length may be inserted
	  *  between these returned left and right sites using a MolModeler object targeting
	  *  the returned IntrctMol.
	  *
	  */
	 IntrctMol mergechains(int firstchain, int secondchain,
			 SiteIdx &leftlinkersite, SiteIdx &rightlinkersite) const;

	 /**
	  * making one chain into two or chop off a terminal segment from one end of a chain
	  *  if newnterm >newcterm+1; the section between the two new terminus will  be discarded.
	  *
	  * if newncterm==-1 the nterm segments is chopped off.
	  *  if newcterm==chainlength, the cterm segments is chooped off
	  *
	  */
	 IntrctMol splitchain(int chainid, int newcterm, int newnterm) const;

	 /**
	  * Build a random peptide segment of a given amino acid sequence and secondary structure sequence and
	  * add that as a new chain of imol_
	  */
	 void buildnewchain(NSPgeometry::XYZ r0, NSPgeometry::XYZ direction,
			 const std::vector<char> &ssseq,
			 const std::vector<std::string> & aaseq=std::vector<std::string>() );
		/**
		 * build a secondary_element with specify direction add by miaoyy
		 */
	 void buildssele(std::vector<char> sseq, std::vector<std::vector<double>> sscrd,std::vector<int> sslen,int loop_option,std::vector<char> &finalss,std::vector<std::string> &ss_info, std::vector<std::string>& ssEH, std::vector<std::string>& ssE, std::vector<std::string>& ssH, std::vector<std::string>& ssC,std::vector<char> direction);

	 /**
     * mv a chain to a new geometric center
     */
	 void movechaincenter(const NSPgeometry::XYZ & newcenter, int chainid=0){
		 target_->translatechaincenter(newcenter*A2NM,chainid);
	 }
	 /**
	  * replace a specified region by a new loop. The coordinates of  input loop are not moved.
	  */
	 void replaceloop
		 	 (const SiteIdx &siteleft, const SiteIdx &siteright,
		 			 const std::vector<IntrctBlck> &newloop){
		int expectedlength=siteright.second-siteleft.second-1;
		assert(newloop.size()==expectedlength);
		assert(newloop[0].getidx2d()==SiteIdx({siteleft.first,siteleft.second+1}));
		for(auto &r:newloop){ target_->pasteblck(r);}
		target_->natoms(); //setup the atom index offsets of IntrctBlcks
		}
		/**
		** output a specified region by SiteIdx.
		**/
		  std::vector<IntrctBlck> outputloop(const SiteIdx &siteleft, const SiteIdx &siteright){
			int expectedlength=siteright.second-siteleft.second-1;
			std::vector<IntrctBlck> newloop; //(expectedlength);
			for (int i = 0; i < expectedlength; i++)
			newloop.push_back(target_->getblck(SiteIdx({siteleft.first,siteleft.second+i+1})));
			return newloop;
		  }
	/**
	 * modify loop conformation
	 *
	 *
	 * The backbone conformation is  generated by LoopConfMaker, which wraps LoopReservoir
	 * if checkclashes_, clashes of the new backbone with the
	 * other part of the target will be checked. The first loop supplied by loopreservoir with
	 * no more than acceptedclashes will be used, or the loop with the least number of clashes
	 * after maxtestnumber tries will be used.
	 * the loop residue length is determined by the loopresiduetype vector.
	 * The new loop residues must be either GLY or ALA, as there is little
	 * point in building other types of side chains on an unoptimized loop.
	 */
//    std::vector<IntrctBlck> newloop(LoopConfMaker &confmaker,
//    		const std::vector<std::string> &loopresiudetypes,
//    		int acceptedclashes=0, int maxtestnumber=100);
    std::vector<IntrctBlck> newloop(LoopConfMaker &confmaker,
    		const std::vector<std::string> &loopresiudetypes,
    		int acceptedclashes=0, int maxtestnumber=100, bool helixinloop=false);
	/**
	 * Determines and returns nearby sites based on CA-CA distance and the distance cutoff.
	 * If the nearby sites have been determined before,
	 * simply returns them.
	 */
	const std::vector<SiteIdx> & getnearbysites(const SiteIdx & site);

	/**
	 * Returns previously determined nearby sites.
	 */
	const std::vector<SiteIdx> & getnearbysites(const SiteIdx & site) const{
		assert(siteneighbors_.find(site)!=siteneighbors_.end());
		return siteneighbors_.at(site);
	}

	/**
	 * Returns previously determined backbone phi-psi at a site.
	 */
	NSPdesignseq::Phipsi getphipsi(const SiteIdx &site) const;

	/**
	 * Calculates and returns the backbone phi-psi at a site.
	 * If pre-calculated values available, simply returns them.
	 */
	NSPdesignseq::Phipsi getphipsi(const SiteIdx &site);

	/**
	 * Returns the local coordinate frame determined by the backbone atoms at a site.
	 *
	 * The returned object can be used to transform sidechain coordinates between local and global
	 * coordinate system.
	 */
	const NSPgeometry::LocalFrame & getlocalframe(const SiteIdx &site) const;

	/**
		 * Calculates the local coordinate frame determined by the backbone atoms at a site.
		 *
		 * If pre-calculated values available, simply returns them.
		 */
	const NSPgeometry::LocalFrame & getlocalframe(const SiteIdx &site);

	/**
	 * returns number of atomic clashes between two residues(IntrctBlck objects).
	 *
	 * If the two blocks are next to each other in sequence (according to the
	 * chain and position indices), covalent exclusion is considered.
	 * The first element of returned array is the number of clashes between mainchain atoms.
	 * The second element is the number of clashes involving sidechain atoms.
	 */
	std::array<int,2> nclashes(const IntrctBlck &blk1,
			const IntrctBlck &blk2) const;
	/**
	 * construct a residue (an IntrctBlck object) of the given residue type.
	 * If multiple rotamers exists for the residue in the rotamer library,
	 * the rotamer that leads to the least number of atomic clashes (if checkclash_)
	 * and the lowest rotamer energy is selected.
	 */
	IntrctBlck buildsidechain(const SiteIdx &site,
				const std::string &sctype){
		return  buildsidechain(site,sctype,checkclash_);
	}
	/**
	 * construct a residue (an IntrctBlck object) of the given residue type.
	 * If multiper rotamers exists for the residue in the rotamer library,
	 * the rotamer that leads to the least number of atomic clashes (if checkclash)
	 * and the lowest rotamer energy is selected.
	 */
	IntrctBlck buildsidechain(const SiteIdx &site,
			const std::string &sctype,bool checkclash);



	/**
	 * make a residue (an IntrctBlck object) use one given rotamer.
	 *
	 * The Block is "floating", i.e., not inserted into the IntrctMol object
	 * pointed to by target_.
	 * this method is called by buildsidechain, which selects between different rotamers.
	 */
	IntrctBlck mkfloatingblck(const SiteIdx & site,NSPdesignseq::Rotamer *rt);

    /**
     * make a floating IntrctBlck object of glysine or alanine residue using the
     * the new coordinates of backbone atoms in backbone site.
     */
    IntrctBlck mkfloatingblck(const SiteIdx &site,
    		const std::string restype, const NSPproteinrep::BackBoneSite &bs) const;

	LoopConfMaker getloopconfmaker(const SiteIdx &leftsite, const SiteIdx &rightsite){
		LoopConfMaker res;
		res.reservoir= mkloopreservoir(leftsite,rightsite);
		res.leftsite=leftsite;
		res.rightsite=rightsite;
		return res;
	}

	/**
	 * returns a loop reservoir that can build closed loops between two sites.
	 * The loop would be at the n terminus if leftsite residue index is negative,
	 * or at the C-terminus if the rightsite residue index >= the tatget chain length
	 */
	LoopReservoir mkloopreservoir(const SiteIdx &leftsite, const SiteIdx &rightsite);
private:
	double neighbor_cutoff2_{1.21}; //!<CA-CA distance cutoff squared,for finding nearby sites
	double atomic_clash_nhb_{0.0841}; //!<atomic clash cutoff squared, for non-hydrogen bonding atom pairs
	double atomic_clash_hb_{0.0576}; //!<atomic clash cutoff squared for hydrogen bonding atom ppairs.
	IntrctMol *target_ { nullptr }; //!<points to the target molecule for modification
	bool checkclash_ { false}; //!<whether to consider atom clashes for rotamer selection and so on
	std::map<SiteIdx,std::vector<SiteIdx>> siteneighbors_; //!<lists of already-computed nearby sites for each site
	std::map<SiteIdx,NSPdesignseq::Phipsi> phipsis_; //!<already-computed phi-psi values
	std::map<SiteIdx, NSPgeometry::LocalFrame> sitelocalframes_; //!<already-computed local coordinate frames
	std::vector<std::vector<NSPproteinrep::BackBoneSite>> backbonesites_;
};
void ssregions(const std::vector<char> &ssseq,std::vector<std::pair<int,int>> &helixregions,
	std::vector<std::pair<int,int>> &strandregions);
}

#endif /* IBLOCK_MOLMODELER_H_ */

/*
 * isolatedpocket.h
 *
 *  Created on: 2020年7月13日
 *      Author: hyiu
 */

#ifndef PROTEINREP_ISOLATEDPOCKET_H_
#define PROTEINREP_ISOLATEDPOCKET_H_
#include "proteinrep/aaconformer.h"
#include "dataio/inputlines.h"
#include <iostream>
#include <memory>
namespace NSPproteinrep{

/**
 * A cluster of amino acid residues and possibly one or more ligand molecules
 *
 * It models a small molecule binding site, such as an enzyme active site.
 */
class IsolatedPocket{
public:
	/**
	 * Uniquely identifies a component of the pocket.
	 *
	 * A component is an AAConfomer object.
	 * It can correspond to a standard natrual amino acid residue or a small molecule.
	 * ComponentID is chain id + integer_residueID +residuetype name
	 * In case of a natural residue component, the residuetype is the capital three letter code.
	 */
	struct ComponentID{
		char chainid{'A'}; //*< different components can be of the same chain id.
		int residueid{0}; //*< different conponents of the same chainid must have different residueid.
		std::string resname{""}; //*< capital three letter code for standard amino acid residues
		ComponentID(){;};
		ComponentID(char cid, int rid, const std::string &rname):
			chainid(cid),residueid(rid),resname(rname){;}
		ComponentID(const AAConformer &aa){
			chainid=aa.chainid_or;
			residueid=aa.residueid_or;
			resname=aa.residuename;
		}
		/**
		 * A helper method.
		 * @param residue_specifier composed of three white space-separated words respectively corresponding to chainid, residue id and residuetype
		 */
		ComponentID(const std::string & residue_specifier){
			std::vector<std::string> words=NSPdataio::parseline(residue_specifier,std::vector<int>());
			chainid=words[0][0];
			residueid=std::stoi(words[1]);
			resname=words[2];
		}

		bool operator <(const ComponentID &b) const{
			if(this->chainid == b.chainid){
				if(this->residueid==b.residueid)
					return this->resname <b.resname;
				else
					return this->residueid <b.residueid;
			} else  return this->chainid <b.chainid;
		}
		bool operator==(const ComponentID &b) const {
			return (!(*this<b) && !(b<*this));
		}
	};
	/**
	 * Uniquely identifies an atom using ComponentID and atom name.
	 *
	 * Atoms are NOT identified by integer indices so that AtomIDs will not be affected whtn
	 *  components are dynamically added or removed.
	 *
	 */
	struct AtomID{
		ComponentID cmpntID;
		std::string atmnm;
		bool operator <(const AtomID &b) const{
				if(this->cmpntID == b.cmpntID){
					return atmnm<b.atmnm;
				} else  return this->cmpntID<b.cmpntID;
			}
		bool operator==(const AtomID &b) const {
				return (!(*this<b) && !(b<*this));
			}
	};
	/**
	 * write the components to os in PDB format
	 */
	std::ostream & writepdb(std::ostream &os);
	/**
	 * extract components from pdb.
	 *
	 * @param pdb the pdb reader that contained all the previously read pdbrecords
	 * @param cmpntids ComponentIDs of residues to extract
	 *                                       If empty or not given, all residues in pdb will be included.
	 */
	void addcomponents(const PdbReader &pdb,const std::vector<ComponentID> & cmpntids=std::vector<ComponentID>());
/**
 * Add an aaconformer as a new component to the isolatedpocket.
 *
 * The new conformer will be given a new chainID and a new residueID.
 * @param newchainid specify the newchainid.
 * The new residueID is automatically determined so that the new conformer' componentID does not overlap with existing ones.
 * @param anchoratoms1
 * @param anchoratomnms
 * @param intcrds
 *  together specifies how the atomic cooridnates of the new component will be translated and rotated.
 *  if anchoratms is not empty, the coordinates of the new conformer will be transformed so that
 * the relative geometries between the  set of three anchoratoms (specified by anchoratomnms)
 *  in the new component relative to the set of three existing anchor atoms (specified by anchoratoms1)
 * are the same as specified by intcrds.
 */
	ComponentID addcomponent(const AAConformer &aaconf, char newchainid,
			const std::vector<AtomID> & anchoratoms1,
			const std::vector<std::string> & anchoratomnms,
			const std::vector<double> &intcrds);
	/**
	 * return a newresidue id that exceeds those already used residue ids for a chain
	 */
	int newresidueid(char newchainid) const {
		int maxusedid=-1;
		for(auto &c:components_) {
			if(c.first.chainid!=newchainid) continue;
			if(c.first.residueid>maxusedid) maxusedid=c.first.residueid;
		}
		return maxusedid<0?1:maxusedid+1;
	}
	bool empty() const{ return components_.empty();}
	std::map<ComponentID,AAConformer> & components(){return components_;}
	const std::map<ComponentID,AAConformer> & components() const {return components_;}
	std::set<ComponentID> &ligandcomponents() {return ligandcomponents_;}
	const std::set<ComponentID> &ligandcomponents() const {return ligandcomponents_;}
private:
	std::map<ComponentID,AAConformer> components_; //*< residues and small molecules comprising the pocket
	std::set<ComponentID> ligandcomponents_;
};
}
#endif /* PROTEINREP_ISOLATEDPOCKET_H_ */

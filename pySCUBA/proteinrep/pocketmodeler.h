/*
 * pocketmodeler.h
 *
 *  Created on: 2020年7月17日
 *      Author: hyiu
 */

#ifndef PROTEINREP_POCKETMODELER_H_
#define PROTEINREP_POCKETMODELER_H_
#include "proteinrep/isolatedpocket.h"
#include "geometry/sphere.h"
namespace NSPproteinrep{
/**
 * For creating, editing and changing conformation of an IsolatedPocket object
 */
class PocketModeler{
public:
	typedef IsolatedPocket::AtomID AtomID;
	typedef IsolatedPocket::ComponentID ComponentID;
	/**
	 * defines a joint or pivot connection that can be used to change the relative
	 * geometries between two parts of the pocket in 6 degrees of freedom.
	 */
	struct Joint{
		bool flexible{true}; //*<if not lexible, the joint will not be used for changing conformations
		std::set<ComponentID> mvcmpnts; //*<specify the components moved by moving around the joint.
		std::vector<AtomID> fx3atms; //*<three atoms (must be non-colinear) from the fixed parts
		std::vector<AtomID> mv3atms; //*<three atoms from the moving parts
		/**
		 * variation ranges (min and max values) of the six internal coordinates defined by the fixed and moving atoms
		 * arranged in the following order (angles and torsions are in radius)
		 * torsion fx3atms[0]-fx3atms[1]-fx3atms[2]-mv3atms[0],
		 * angle fx3atms[1]-fx3atms[2]-mv3atms[0]
		 * distance fx3atms[2]-mv3atms[0]
		 * angle fx3atms[2]-mv3atms[0]-mv3atms[1]
		 * torsion fx3atms[1]-fx3atms[2]-mv3atms[0]-mv3atms[1]
		 * torsion fxatms[2]-mv3atms[0]-mvatms[1]-mv3stms[2]
		 */
		std::vector<std::pair<double,double>> intcrdrngs;
	};
	/**
	 * defines natural amino acid residue component whose rotamer state can be varied to
	 * change the conformation of the pocket.
	 * Three given atoms (not necessarily backbone atoms) from the residue will be fixed with respect
	 * the rest of the pocket in conformation change
	 */
	struct FlexRotamer{
		ComponentID cmpntid;
		std::vector<std::string> fxatoms; //*< three fixed non-colinear atoms of the residue
	};
	/**
	 * create a pocket in a process controlled interactively.
	 *
	 * This function will
	 * allocate the pocket_
	 * read  components from PDB
	 * add amino amino acid residues
	 * define the geometries of the added residue relative to the other components by defining Joint and
	 * FlexRotamer
	 */
	void interactivesetup(std::ostream &os, std::istream &is);
/*	void addjoint(bool isflexible,const std::vector<AtomID> &fx3atms,
			const std::vector<AtomID> &mv3atms,
			const std::vector<std::pair<double,double>> & intcrdrngs,
			std::set<ComponentID> mvcmpnts);*/

	/**
	 * change the pocket conformation by a random move around a joint
	 */
	void mvaroundjoint(const Joint & joint);

	/**
	 * change the pocket conformation by random move around a joint
	 * trial random moves are repeated until the number of atomic clashes are
	 * no more than the starting number of clashes (given by nclashes)
	 * @param nclashes the starting number of clashes. Delivered with the new number of clashes
	 */
	void mvaroundjoint(const Joint &joint, int &nclashes){
		int ntotal=nclashes+1;
		while(ntotal>nclashes){
			mvaroundjoint(joint);
			ntotal=totalclashes();
		}
		nclashes=ntotal;
	}
	/**
	 * change the pocket conformation by change a flexible rotamer
	 *
	 */
	void mvrotamer(const FlexRotamer & flexrt);

	/**
	 * change the pocket conformation by random move around a joint
	 * trial random moves are repeated until the number of atomic clashes are
	 * no more than the starting number of clashes (given by nclashes)
	 * @param nclashes the starting number of clashes. Delivered with the new number of clashes
	 */
	void mvrotamer(const FlexRotamer &flexrt, int &nclashes){
		int ntotal=nclashes+1;
		while(ntotal>nclashes){
			mvrotamer(flexrt);
			ntotal=totalclashes();
		}
		nclashes=ntotal;
	}
    /**
     * generate a clash-free trial conformation
     *
     * It first considers random moves of all joints and flexrotamers,
     * and then randomly change individual joints and
     * flexrotamers to reduce the number of clashes until a clash-free
     * conformation is obtained.
     * If the number of individual moves exceeded a certain number without
     * fielding a clash-free conformation, it will return false.
     */
	bool genconformation_trial();

	/**
	 * call genconformation_trial repeatedly to obtain a clash-free conformation
	 *
	 * @param maxtrial specifies the maximum number of calling genconformation_trial.
	 */
	bool genconformation(int maxtrial){
		int ntrial=0;
		bool success=false;
		while(!success && ntrial++ <maxtrial){
			success=genconformation_trial();
		}
		return success;
	}
	NSPgeometry::XYZ ligandcenter() const;
	void rigidtranslation(const NSPgeometry::XYZ & shift);
	void rigidmv(double maxtrans,double maxrotate, const NSPgeometry::Sphere &ligandcentersphere);
	bool randominternalmv(int &nclashes);
	std::vector<NSPgeometry::XYZ> copycrdto() const;
	void copycrdfrom(const std::vector<NSPgeometry::XYZ> &crds);
	void writepdb(std::ostream &os){pocket_->writepdb(os);}
	IsolatedPocket &pocket() {return *pocket_;}
	const IsolatedPocket &pocket() const {return *pocket_;}
	std::vector<Joint> & joints(){return joints_;}
	const std::vector<Joint> &joints() const {return joints_;}
	std::vector<FlexRotamer> & flexrotamers(){return flexrotamers_;}
	const std::vector<FlexRotamer> & flexrotamers() const {return flexrotamers_;}
	int totalclashes();
private:
	std::shared_ptr<IsolatedPocket> pocket_;
	std::vector<Joint> joints_;
	std::vector<FlexRotamer> flexrotamers_;
	std::map<AtomID,std::set<AtomID>> clashexclusions_; //atom pairs excluded from clash checking
	std::map<ComponentID,std::set<ComponentID>> cexclusions_;//component pairs excluded from clash checking
	bool chainup2joints(Joint & j1, Joint &j2);
	void chainupjoints();
	void setclashexclusions();
	int nclashes(const ComponentID &c1, const ComponentID &c2);

};
}



#endif /* PROTEINREP_POCKETMODELER_H_ */

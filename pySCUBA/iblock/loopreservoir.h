/*
 * loopresevoir.h
 *
 *  Created on: 2019年12月21日
 *      Author: hyiu
 */

#ifndef IBLOCK_LOOPRESERVOIR_H_
#define IBLOCK_LOOPRESERVOIR_H_
#include "backbone/backbonesite.h"
#include "dstl/randomengine.h"
namespace NSPintrct{
/**
 *constructs and stores a reservoir of closed loop conformations
 *
 *The N-terminal backbone site and the C-terminal backbone site are all fixed.
 *The N and CA atoms of the first loop residue are at fixed positions
 *because of the fixed CA-C=O atoms of the previous residue and the rigid backbone
 *geometry. The CA-C=O atoms of the last loop residue is also at fixed positions.
 *
 *If the loop is at the N-terminus or C-terminus, there will be no need to close.
 */
class LoopReservoir{
public:
	typedef std::vector<NSPproteinrep::BackBoneSite> LoopConf; //*<loop conformation
	LoopReservoir(){;}
	LoopReservoir(const NSPproteinrep::BackBoneSite &leftend,
				const NSPproteinrep::BackBoneSite &rightend,int randomseed=0){
		init(leftend,rightend,randomseed);
	}
	/**
	 * Just initialize the fixed terminus.
	 * No loop built at construction. The current implementation
	 *
	 * @para leftend the N-terminal backbone site
	 * @para rightend the C-terminal backbone site
	 * @para randomseed if randomseed !=0, the global randomengine instance will be reseeded
	 * WARNING: this version calls backbonebuilder which
	 * makes use of the GLOBAL instance of the randomengine.
	 * So if randomseed !=0, the GLOBAL random engine instance will be reseeded.
	 */
	void init(const NSPproteinrep::BackBoneSite &leftend,
			const NSPproteinrep::BackBoneSite &rightend,int randomseed=0)
	{
		static int defalutseed=1357;
		anchor_left_=leftend;
		anchor_right_=rightend;
		if(randomseed!=0)
				NSPdstl::RandomEngine<>::getinstance().reseed(randomseed);
	}

	void init_nterminus(const NSPproteinrep::BackBoneSite &rightend,
			int randomseed=0){
		nterminal_=true;
		anchor_right_=rightend;
		if(randomseed!=0)
				NSPdstl::RandomEngine<>::getinstance().reseed(randomseed);
	}
	void init_cterminus(const NSPproteinrep::BackBoneSite &leftend,
			int randomseed=0){
		cterminal_=true;
		anchor_left_=leftend;
		if(randomseed!=0)
				NSPdstl::RandomEngine<>::getinstance().reseed(randomseed);
	}
	/**
	 * return the next closed loop conformation for the given length
	 *
	 * if no next conformation is available, will try to generate a collection of
	 * closed loops and put the new solutions in the reservior.
	 * If the closing try failed, an empty vector will be returned.
	 */
	std::shared_ptr<LoopConf> popconf(int length);
	// especially used for initially establish helix_PhiPsi within loops.
	std::shared_ptr<LoopConf> popconf(int length, bool hashelix);

private:
	/**
	 * stores preconstructed loop conformations of different lenghths
	 */
	std::map<int,std::vector<std::shared_ptr<LoopConf>>> reservoir_;
	int max_tryclosing_{100};
	bool nterminal_{false}; //<* n-terminal loop
	bool cterminal_{false}; //<* c-terminal loop
	/**
	 * for each loop length, points to the location of a new loop to return
	 * on next popconf call.
	 */
	std::map<int,int> posi_;
	NSPproteinrep::BackBoneSite anchor_left_;//<* fixed N-terminal site
	NSPproteinrep::BackBoneSite anchor_right_;//<* fixed C-terminal site
//todo	NSPdstl::RandomEngine<> rng_;
};

}

#endif /* IBLOCK_LOOPRESERVOIR_H_ */

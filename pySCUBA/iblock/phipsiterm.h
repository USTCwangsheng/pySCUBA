/*
 * phipsiterm.h
 *
 *  Created on: 2019年12月4日
 *      Author: hyiu
 */

#ifndef IBLOCK_PHIPSITERM_H_
#define IBLOCK_PHIPSITERM_H_
#include "pdbstatistics/phipsidistr.h"
#include "iblock/phipsicodes.h"

namespace NSPintrct{
struct IntrctPara;
/**
 * wraps two different ways of calculating the phipsi energy
 *
 * For teminal residues, direct interpolation of the 1-D statistical distribution
 * FOr middle residues, neuronetwork representing the 2-D energy surface
 */
struct PhiPsi {
	enum MODE {USUAL,NTERM,CTERM};
	const NSPpdbstatistics::PhiPsiDistr *distr;
//	std::vector<int> atoms;
	int mode;
	int siteid{-1};
	int chainid{0};
	PhiPsi(int sid,int cid,bool isnterm=false,bool iscterm=false):siteid(sid),chainid(cid){
		if(isnterm) mode=NTERM;
		else if(iscterm) mode=CTERM;
		else mode=USUAL;
		distr=&(NSPpdbstatistics::PhiPsiDistr::coildistr());
	}
	PhiPsi(int sid,int cid,std::string resname, std::string resnext,bool isnterm=false,bool iscterm=false):siteid(sid),chainid(cid){
		if(isnterm) mode=NTERM;
		else if(iscterm) mode=CTERM;
		else mode=USUAL;
		distr=&(NSPpdbstatistics::PhiPsiDistr::phipsidistr(resname,resnext));
	}
	double energy(const IntrctPara &para,const std::vector<PhiPsiCodes> &phipsicodes,
			std::vector<NSPgeometry::XYZ> *forces) const;
};

}


#endif /* IBLOCK_PHIPSITERM_H_ */

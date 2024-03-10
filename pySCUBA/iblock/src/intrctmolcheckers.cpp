/*
 * intrctmolcheckers.cpp
 *
 *  Created on: 2019年12月26日
 *      Author: hyiu
 */
#include "iblock/intrctmol.h"

using namespace NSPintrct;

bool IntrctMol::checkblckindices() const{
	for(int i=0;i<nchains();i++){
		for(int j=0;j<nresidues(i);j++){
			if( (*this)({i,j}).getidx2d() != NSPdstl::Idx2D({i,j})){
				return false;
			}
		}
	}
	return true;
}
bool IntrctMol::checkaoffsets() const {
	int aoff=0;
	for(int i=0;i<nchains();i++){
		for(int j=0;j<nresidues(i);j++){
			if((*this)({i,j}).aoffset !=aoff) return false;
			aoff += (*this)({i,j}).natoms();
		}
	}
	return true;
}

bool IntrctMol::checkbsinchains() const {
	for(auto &bss:bsinchains){
		for (auto & bs:bss){
			auto & blk=(*this)({bs.chainid,bs.resid});
			if(blk.atomname(bs.nid-blk.aoffset) != "N") return false;
			if(blk.atomname(bs.caid-blk.aoffset) != "CA") return false;
			if(blk.atomname(bs.cid-blk.aoffset) != "C") return false;
			if(blk.atomname(bs.oid-blk.aoffset) != "O") return false;
		}
	}
	return true;
}

bool IntrctMol::checkscinchains() const {
	for(int c=0;c<bsinchains.size();c++){
		for(int r=0;r<bsinchains[c].size();r++){
			auto &blk=(*this)({bsinchains[c][r].chainid,bsinchains[c][r].resid});
			if(scinchains[c][r].restype != blk.resname()) return false;
		}
	}
	return true;
}

bool IntrctMol::checkcrds() const {
	for(int i=0;i<nchains();++i){
		for(int j=0;j<nresidues(i);j++){
			auto &blk=(*this)({i,j});
			for(int a=0;a<blk.natoms();a++){
				if((crds_all_[a+blk.aoffset]-blk.crds[a]).squarednorm() > 1.e-8) return false;
			}
		}
	}
	return true;
}

bool IntrctMol::checkatomsfixed() const {
	for(int i=0;i<nchains();++i){
		for(int j=0;j<nresidues(i);j++){
			auto &blk=(*this)({i,j});
			for(int a=0;a<blk.natoms();a++){
				bool fixed=false;
				if(blk.activemod==IntrctBlck::INACTIVE) fixed=true;
				else if(blk.activemod==IntrctBlck::SIDECHAIN)
					if(blk.topo->atomips[a].ismainchain) fixed=true;
				if(isatomfixed_[a+blk.aoffset]!=fixed) return false;
			}
		}
	}
	return true;
}

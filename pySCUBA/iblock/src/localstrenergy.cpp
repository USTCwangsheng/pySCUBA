/*
 * localstrenergy.cpp
 *
 *  Created on: 2019年12月4日
 *      Author: hyiu
 */
#include "iblock/intrctmol.h"
#include "iblock/nnterms.h"
using namespace NSPintrct;
void IntrctMol::localstrenergy(const IntrctPara &param,
		const std::vector<BSInChain> & bsinchain,
		const std::vector<PhiPsiCodes> & phipsicodes,
		std::vector<NSPgeometry::XYZ> *dedx_sys) const {
	double w=param.weight_localstr;
	for (int i = 2; i < phipsicodes.size() - 2; ++i) {
		if(!windowactive(bsinchain,i,2)) continue;
		LSNNTerm lsnnterm;
		std::vector<DvDxi> dedx;
		lsnnterm.setup(phipsicodes, i);
		double e = w*lsnnterm.outvalue(&dedx);
		if(param.enedetails){
			Results::ofstreams(param.jobname,"LocalStrEne")
			   <<"chain "<<bsinchain[i].chainid
			   <<"residue "<<i <<" e = "<<e<<std::endl;
		}
		mol_(bsinchain[i].chainid,bsinchain[i].resid).energies[IntrctBlck::LOCALSTR]=e;
		for (auto &d : dedx) {
			(*dedx_sys)[d.first] = (*dedx_sys)[d.first] + w * d.second;
		}
	}
}
void IntrctMol::forces_localstr(const IntrctPara &param,std::vector<NSPgeometry::XYZ> &dedx) const{
	if(phipsicodes.empty()) calcphipsicodes();
	auto &chains = mol_.D_;
	for (auto &c : chains) {
		for (auto &r : c)
			r.energies[IntrctBlck::LOCALSTR] = 0.0;
	}
	for(int i=0;i<bsinchains.size();++i){
		localstrenergy(param,bsinchains[i],phipsicodes[i],&dedx);
	}
}



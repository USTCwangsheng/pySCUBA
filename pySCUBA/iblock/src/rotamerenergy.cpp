/*
 * rotamerenergy.cpp
 *
 *  Created on: 2019年12月11日
 *      Author: hyiu
 */
#include "iblock/intrctmol.h"
#include "iblock/nnterms.h"
using namespace NSPintrct;
void IntrctMol::forces_rotamer(const IntrctPara &param, std::vector<NSPgeometry::XYZ> &dedx_sys)const {
	if(conformercodes.empty())calcconformercodes();
	for (auto &c : mol_.D_) {
		for (auto & r : c)
			r.energies[IntrctBlck::ROTAMER] = 0.0;
	}
	int nchains = phipsicodes.size();
	double ene = 0.0;
	int pidx = -1;
	for (int c = 0; c < nchains; ++c) {
		for (int p = 0; p < phipsicodes[c].size(); ++p) {
			++pidx;
			int resid=bsinchains[c][p].resid;
			const IntrctBlck &blk=mol_(c,resid);
			if(blk.activemod>IntrctBlck::SIDECHAIN){
				if(!windowactive(bsinchains[c],p,1)) continue;
			}
			if (scinchains[c][p].kaiatoms.empty())
				continue;
			std::vector<DvDxi> dedx;
			double e;
			if(p==0 || p==conformercodes[c].size()-1){
				NN_KaiTerm_T kaiterm(conformercodes[c][p]);
				e=kaiterm.outvalue(&dedx);
			} else {
				NN_KaiTerm kaiterm(conformercodes[c][p]);
				e = kaiterm.outvalue(&dedx);
			}
				double w=param.weight_rotamer;
			blk.energies[IntrctBlck::ROTAMER] = w*e;
			if(param.enedetails){
				Results::ofstreams(param.jobname,"rotamerEne")
				   <<"chain "<<blk.getidx2d().first
				   <<"residue "<<blk.getidx2d().second
				   <<" e = "<<e <<std::endl;
			}
			for (auto &d : dedx) {
				dedx_sys.at(d.first) = dedx_sys.at(d.first) + w * d.second;
			}
		}
	}
}




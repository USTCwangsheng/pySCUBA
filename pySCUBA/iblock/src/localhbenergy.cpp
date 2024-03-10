/*
 * localhbenergy.cpp
 *
 *  Created on: 2019年12月12日
 *      Author: hyiu
 */
#include "iblock/intrctmol.h"
#include "iblock/nnterms.h"
using namespace NSPintrct;
void IntrctMol::forces_localhb(const IntrctPara &param,
		std::vector<NSPgeometry::XYZ> & dedx_sys) const {
	assert(offsetsok_);
	auto &chains = mol_.D_;
	for (auto &c : chains) {
		for (auto &r : c)
			r.energies[IntrctBlck::LOCALHB] = 0.0;
	}
	if (param.weight_localhb == 0) return;
	int minsep = 5;
	double ene = 0.0;
	int chainid = -1;
	for (auto &bss : bsinchains) {
		chainid++;
		for (int s1 = 1; s1 < bss.size() - 3; s1++) {
			auto & bs1 = bss[s1];
			if(!windowactive(bss,s1+2,3)) continue;
			auto &blk1 = mol_(bs1.chainid, bs1.resid);
			for (int s2 = s1 + 3; s2 <= s1 + minsep; ++s2) {
				if (s2 >= bss.size())
					continue;
				auto &bs2 = bss[s2];
				auto &blk2 = mol_(bs2.chainid, bs2.resid);
				double w = blk1.weight_mc* blk2.weight_mc * param.weight_localhb;
				if (w == 0) continue;
				std::vector<int> atoms;
				atoms.push_back(bss[s1 - 1].cid);
				atoms.push_back(bss[s1 - 1].oid);
				atoms.push_back(bss[s1].nid);
				atoms.push_back(bss[s2 - 1].cid);
				atoms.push_back(bss[s2 - 1].oid);
				atoms.push_back(bss[s2].nid);
				std::vector<DvDxi> dvdx;
				LocalBbHBGeoNNTerm nnmodel;
				nnmodel.setup(crds_all_, atoms);
				double ehb = w * nnmodel.outvalue(&dvdx);
				if (ehb == 0)
					continue;
				if (ehb != 0.0) {
					for (auto & d : dvdx)
						dedx_sys.at(d.first) = dedx_sys.at(d.first)
								+ w * d.second;
				}
				if (param.enedetails) {
					Results::ofstreams(param.jobname, "LocalHBEne") << "chain "
							<< chainid << "residues " << s1 << "-" << s2
							<< " e = " << ehb << std::endl;
				}
				blk1.energies[IntrctBlck::LOCALHB]+= 0.5 * ehb;
				blk2.energies[IntrctBlck::LOCALHB]+= 0.5 * ehb;
			}
		}
	}
}


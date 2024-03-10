/*
 * sitepaitenergy.cpp
 *
 *  Created on: 2019年12月8日
 *      Author: hyiu
 */
#include "iblock/intrctmol.h"
#include "iblock/nnterms.h"
using namespace NSPintrct;
using namespace NSPgeometry;
bool sitepairactive(const IntrctBlck &blk1, const IntrctBlck &blk2) {
	static const int minsep = 5;
	static const double rcut2=0.9025;
	NSPdstl::Idx2D i2d1 = blk1.getidx2d();
	NSPdstl::Idx2D i2d2 = blk2.getidx2d();
	if (i2d1.first == i2d2.first) {
		int dp = i2d2.second - i2d1.second;
		if (dp <= minsep && dp >= -minsep)
			return false;
	}
	bool a1 = blk1.getmol()->mainchainwindowactive(i2d1, 1);
	bool a2 = blk2.getmol()->mainchainwindowactive(i2d2, 1);
	if(a1|| a2) {
		double r2=(blk1.crds[1]-blk2.crds[1]).squarednorm();
		return r2<rcut2;
	}
	return false;
}
void IntrctMol::forces_sitepair(const IntrctPara &param,
		std::vector<NSPgeometry::XYZ> &dedx_sys) const {
	if(phipsicodes.empty()) calcphipsicodes();
	for (auto &c : mol_.D_) {
		for (auto & r : c)
			r.energies[IntrctBlck::SITEPAIR] = 0.0;
	}
	int nchains = bsinchains.size();
	for (int n = 0; n < nchains; ++n) {
		const std::vector<BSInChain> &bsn = bsinchains[n];
		for (int m = n; m < nchains; ++m) {
			const std::vector<BSInChain> &bsm = bsinchains[m];
			for (int i = 1; i < bsn.size() - 1; ++i) {
				int jstart = 1;
				int jend = bsm.size() - 1;
				if (m == n) {
					jstart = i + 1;
				}
				const IntrctBlck & blk1 = mol_(bsn[i].chainid, bsn[i].resid);
				for (int j = jstart; j < jend; ++j) {
					if (m == n) {
						if (samesssegment(sscodes[n], i, j))
							continue;
					}
					const IntrctBlck &blk2 = mol_(bsm[j].chainid, bsm[j].resid);
					if (!sitepairactive(blk1, blk2))
						continue;
					SitePairNNTerm spterm;
					spterm.setup(crds_all_d_, bsn, phipsicodes[n], i, bsm,
							phipsicodes[m], j);
					double sw = blk1.weight_mc * blk2.weight_mc
							* param.weight_sitepair;
					std::vector<DvDxi> dedx;
					double eh = 0.5 * sw*spterm.outvalue(&dedx);
					blk1.energies[IntrctBlck::SITEPAIR] += eh;
					blk2.energies[IntrctBlck::SITEPAIR] += eh;
					if(param.enedetails){
						Results::ofstreams(param.jobname,"sitepairEne")
						   <<"chain "<<bsn[i].chainid
						   <<"residue "<<bsn[i].resid <<":"
						   <<"chain "<<bsm[j].chainid
						   <<"residue "<<bsm[j].resid
						   <<" e = "<<2.0*eh<<std::endl;
					}
					for (auto &d : dedx) {
						dedx_sys.at(d.first) = dedx_sys.at(d.first)
								+ sw * d.second;
					}
				}
			}
		}
	}
}


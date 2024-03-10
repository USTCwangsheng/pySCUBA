/*
 * shakebonds.cpp
 *
 *  Created on: 2019年12月16日
 *      Author: hyiu
 */

#include "sampling/shakebonds.h"
#include "iblock/intrctmol.h"
#include <cmath>
#include <cassert>
#include <functional>
using namespace NSPsampling;
using namespace NSPintrct;
bool ShakeBonds::shake(const std::vector<double> & masses,
		const std::vector<double> &crdref, std::vector<double> &crd) const {
	if(!ison_) return true;
	int natoms = crd.size() / ndim_;
	std::vector<bool> moved(natoms, true);
	std::vector<bool> moved_next(natoms, false);
	double tolerance = 1.e-4;
	bool done = false;
	while (!done) {
		done = true;
		for (int ib = 0; ib < bonds.size(); ++ib) {
			int i = bonds[ib].i;
			int j = bonds[ib].j;
			if (!moved[i] && !moved[j])
				continue;
			std::vector<double> r = rij(crd, i, j);
			double b2 = norm2(r);
			double diff2 = bonds[ib].b02 - b2;
			if (fabs(diff2) <= tolerance * bonds[ib].b02)
				continue;
			std::vector<double> rref = rij(crdref, i, j);
			double sp = dot(r, rref);
			if (sp < 1.e-10)
				return false;
			double invmassi = 1.0 / (masses[i * ndim_]);
			double invmassj = 1.0 / (masses[j * ndim_]);
			double lamda = diff2 / (2.0 * sp * (invmassi + invmassj));
			for (int d = 0; d < ndim_; ++d) {
				crd[i * ndim_ + d] -= lamda * rref[d] * invmassi;
				crd[j * ndim_ + d] += lamda * rref[d] * invmassj;

			}
			moved_next[i] = true;
			moved_next[j] = true;
			done = false;
		}
		moved = moved_next;
		moved_next.assign(moved_next.size(), false);
	}
	return true;
}

std::shared_ptr<ShakeBonds> NSPsampling::make_shakebonds(IntrctMol &mol){
	std::shared_ptr<ShakeBonds> sbds_p=std::shared_ptr<ShakeBonds>(new ShakeBonds());
	ShakeBonds & sbds=*sbds_p;
	std::vector<bool> atomfixed=mol.isfixedvec();
	for (int ic=0;ic<mol.nchains();++ic){
		for(int ir=0;ir<mol.nresidues(ic);++ir){
			const IntrctBlck &blk=mol.getblck(NSPdstl::Idx2D{ic,ir});
			NSPdstl::Idx2D idx2d=blk.getidx2d();
			for(auto &bnd:blk.topo->bonds){
				bool exist=true;
				for(int i=0;i<2;++i) {
					AtomIdx3D i3d=bnd.aidx.at(i);
					if(!(i3d[0]==0 && i3d[1]==0)){
						NSPdstl::Idx2D bidx2d{idx2d.first+i3d[0],idx2d.second+i3d[1]};
						if(bidx2d.second>=mol.nresidues(bidx2d.first)) {
							exist=false;
						}
					}
				}
				if(!exist) continue;
				int iidx=blk.getatomindex_abs(bnd.aidx[0]);
				int jidx=blk.getatomindex_abs(bnd.aidx[1]);
				if(atomfixed[iidx] && atomfixed[jidx]) continue;
				sbds.addbond(iidx,jidx,bnd.v0);
			}
		}
	}
	return sbds_p;
}


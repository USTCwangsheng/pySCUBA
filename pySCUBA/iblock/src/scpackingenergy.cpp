/*
 * scpackingenergies.cpp
 *
 *  Created on: 2019年12月12日
 *      Author: hyiu
 */
#include "iblock/intrctmol.h"
#include "iblock/vsctypes.h"
using namespace NSPintrct;
void IntrctMol::forces_scpacking(const IntrctPara &param,
		std::vector<NSPgeometry::XYZ> &dedx) const {
	assert(offsetsok_);
	auto &chains = mol_.D_;
	int residue = 0;
	for (auto &c : chains) {
		for (auto &r : c)
			r.energies[IntrctBlck::SCPACKING] = 0.0;
	}
	for (auto &c : chains) {
		for (auto &r : c) {
			r.scpackingforces(param, &dedx);
		}
	}
}
double IntrctBlck::scpackingforces(const IntrctPara &param,
		std::vector<NSPgeometry::XYZ> *dedx_sys) const {
	double esum = 0.0;
	for (auto & entry : nblist_S) {
		double ene = blkscpackingforces(param, (*molsystm_)[entry.first],
				entry.second, dedx_sys);
		this->energies[SCPACKING] += 0.5 * ene;
		(*molsystm_)[entry.first].energies[SCPACKING] += 0.5 * ene;
		esum += ene;
	}
	return esum;
}
double IntrctBlck::blkscpackingforces(const IntrctPara &param,
		const IntrctBlck &blk2, const AtomNbList &nblist,
		std::vector<NSPgeometry::XYZ> *dedx_sys) const {
	double esum = 0.0;
	assert(nblist.size() == natoms());
	double w = param.weight_scpacking;
	for (int i = 0; i < natoms(); ++i) {
		int iidx = aoffset + i;
		bool isoft = softsc && topo->atomips[i].issidechain;
		for (auto j : nblist.at(i)) {
			int jidx = blk2.aoffset + j;
			std::vector<NSPgeometry::XYZ> drdx;
			double r = NSPgeometry::distance(crds[i], blk2.crds[j], &drdx);
			double dedr;
			double ene;
			bool jsoft = blk2.softsc && blk2.topo->atomips[j].issidechain;
			if (isoft || jsoft) {
				ene = w
						* VSCType::softpackingenefuncs.getenefunc(
								topo->atomips[i].sterictype,
								blk2.topo->atomips[j].sterictype).energy(r,
								&dedr);
			} else {
				ene = w
						* VSCType::packingenefuncs.getenefunc(
								topo->atomips[i].sterictype,
								blk2.topo->atomips[j].sterictype).energy(r,
								&dedr);
			}
			if (param.enedetails) {
				Results::ofstreams(param.jobname, "SCPackingEne")
						<< atomstring(i) << "-" << blk2.atomstring(j) << " e = "
						<< ene << std::endl;
			}
			(*dedx_sys)[iidx] = (*dedx_sys)[iidx] + w * dedr * drdx[0];
			(*dedx_sys)[jidx] = (*dedx_sys)[jidx] + w * dedr * drdx[1];
			esum += ene;
		}
	}
	return esum;
}


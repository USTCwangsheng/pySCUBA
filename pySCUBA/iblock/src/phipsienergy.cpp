/*
 * phipsienergy.cpp
 *
 *  Created on: 2019年12月4日
 *      Author: hyiu
 */
#include "iblock/phipsiterm.h"
#include "iblock/intrctparam.h"
#include "iblock/intrctmol.h"
#include "iblock/nnterms.h"
using namespace NSPintrct;
using namespace NSPgeometry;
double PhiPsi::energy(const IntrctPara &para,
		const std::vector<PhiPsiCodes> &phipsicodes,
		std::vector<NSPgeometry::XYZ> *dedx_sys) const {
	static const double rad = 180.0 / 3.14159265;
	std::vector<XYZ> dphidr;
	std::vector<XYZ> dpsidr;
	double phipsiwght = para.weight_phipsi;
	auto & codes = phipsicodes[siteid];
	double dedpsi = 0.0;
	double dedphi = 0.0;
	double ene;
	if (mode == NTERM) {
		ene = distr->intplene_psi(codes.psi * rad, &dedpsi);
		dedpsi *= phipsiwght * rad;
		for (auto &d : codes.dpsidx)
			(*dedx_sys)[d.first] = (*dedx_sys)[d.first] + dedpsi * d.second;
	} else if (mode == CTERM) {
		ene = distr->intplene_phi(codes.phi * rad, &dedphi);
		dedphi *= phipsiwght * rad;
		for (auto &d : codes.dphidx)
			(*dedx_sys)[d.first] = (*dedx_sys)[d.first] + dedphi * d.second;
	} else {
		PhiPsiNNTerm<PhiPsiNNTerm<>::MIXCOIL> phipsinnterm;
		std::vector<DvDxi> dedxi;
		ene = phipsinnterm.phipsiene(phipsicodes, siteid, &dedxi);
		for (auto &d : dedxi) {
			(*dedx_sys)[d.first] = (*dedx_sys)[d.first] + phipsiwght * d.second;
		}
	}
	return ene * phipsiwght;
}
void IntrctMol::forces_phipsi(const IntrctPara &param,
		std::vector<NSPgeometry::XYZ> &dedx) const {
	if(phipsicodes.empty()) calcphipsicodes();
	auto &chains = mol_.D_;
	for (auto &c : chains) {
		for (auto &r : c)
			r.energies[IntrctBlck::PHIPSI] = 0.0;
	}
	int cid = 0;
	for (auto & bss : bsinchains) {
		for (int s = 0; s < bss.size(); ++s) {
			if (!windowactive(bss, s, 1))
				continue; //center and neighboring backbones not active
			PhiPsi pp(s, cid, s == 0, s == bss.size() - 1);
			double ene = pp.energy(param, phipsicodes[cid], &dedx);
			getblck(NSPdstl::Idx2D(bss[s].chainid, bss[s].resid)).energies[IntrctBlck::PHIPSI] =
					ene;
			if (param.enedetails) {
				Results::ofstreams(param.jobname, "PhiPsiEne") << "chain "
						<< cid << " posi " << s << " Phi,Psi = "
						<< phipsicodes[cid][s].phi * 180.0 / 3.14159265 << ","
						<< phipsicodes[cid][s].psi * 180.0 / 3.14159265
						<< "  e = " << ene << std::endl;
			}
		}
		++cid;
	}
}

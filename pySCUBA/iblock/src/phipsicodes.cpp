/*
 * phipsicodes.cpp
 *
 *  Created on: 2019年12月4日
 *      Author: hyiu
 */
#include "iblock/intrctmol.h"
#include "iblock/phipsicodes.h"
#include "designseq/ResName.h"
#include "iblock/nnterms.h"
using namespace NSPintrct;
std::vector<std::vector<BSInChain>> NSPintrct::makebsinchains(
		const IntrctMol &mol) {
	std::vector<std::vector<BSInChain>> bschains(mol.nchains());
//	auto &rn = NSPdesignseq::ResName::resname();
	for (int c = 0; c < mol.nchains(); ++c) {
		auto & bsc = bschains[c];
		for (int r = 0; r < mol.nresidues(c); ++r) {
			std::vector<int> atoms =
					mol.getblck(NSPdstl::Idx2D { c, r }).getbackboneatoms();
			if (atoms.size() == 4) {
				bool active=(mol.getblck(NSPdstl::Idx2D(c,r)).activemod==IntrctBlck::ALL);
				bsc.push_back(
						BSInChain(atoms[0], atoms[1], atoms[2], atoms[3],active));
				bsc.back().chainid=c;
				bsc.back().resid=r;
			}
		}
	}
	return bschains;
}
bool NSPintrct::windowactive(const std::vector<BSInChain> &bsinchain, int center,int halfw){
	int first=center-halfw;
	if(first<0) first=0;
	int last=center+halfw+1;
	if(last>=bsinchain.size()) last=bsinchain.size();
	for(int i=first;i<last;++i) if(bsinchain[i].active) return true;
	return false;
}
bool NSPintrct::samesssegment(const std::vector<SSCode> &sscodes,int p1, int p2){
	assert(p1<p2);
	assert(p1<sscodes.size());
	assert(p2<sscodes.size());
	int ssid1 = sscodes[p1].ssid;
	int ssid2 = sscodes[p2].ssid;
	if (ssid1 == NN_SSTerm::COIL || ssid2 == NN_SSTerm::COIL)
		return false;
	int ssid_p = ssid1;
	if (ssid_p == NN_SSTerm::TERMINUS)
		ssid_p = sscodes[p1 + 1].ssid;
	for (int p = p1 + 1; p < p2; ++p) {
		if (sscodes[p].ssid != ssid_p)
			return false;
	}
	if (ssid2 != ssid_p && ssid2 != NN_SSTerm::TERMINUS)
		return false;
	return true;
}
std::vector<PhiPsiCodes> NSPintrct::makephipsicodes(
		const std::vector<NSPgeometry::XYZ> &crds,
		const std::vector<BSInChain> & bsinchain) {
	int nsite = bsinchain.size();
	std::vector<PhiPsiCodes> result(nsite);
#pragma omp parallel for schedule(dynamic,1)
	for (int i = 0; i < nsite; ++i) {
		//	for(auto it=bsinchain.begin(); it!=bsinchain.end();++it){
		auto it = bsinchain.begin() + i;
		result[i] = PhiPsiCodes(crds, it, it == bsinchain.begin(),
				it == bsinchain.begin() + bsinchain.size() - 1);
	}
	return result;
}
std::vector<SSCode> NSPintrct::estimatess(const std::vector<PhiPsiCodes> &phipsicodes){
	int nsite=phipsicodes.size();
	std::vector<SSCode> result(nsite);
#pragma omp parallel for schedule(dynamic,1)
	for(int i=0;i<nsite;++i) {
		auto & sc=result.at(i);
		sc.p3=NN_SSTerm().probabilities(phipsicodes,i,&(sc.dp3dx));
		sc.ssid=NN_SSTerm::sstype(sc.p3);
	}
	return result;
}

std::vector<double> PhiPsiCodes::gettriangcodes(double ang,
		const std::vector<DvDxi> &dadx, std::vector<std::vector<DvDxi>> *dcdx) {
	const std::vector<double> nang { 1.0, 2.0, 4.0, 8.0 };
	std::vector<double> res;
	dcdx->clear();
	for (double n : nang) {
		double c = cos(n * ang);
		double s = sin(n * ang);
		res.push_back(c);
		dcdx->push_back(std::vector<DvDxi>());
		auto & dccdx = dcdx->back();
		for (auto &d : dadx)
			dccdx.push_back(-(s * n) * d);
		res.push_back(s);
		dcdx->push_back(std::vector<DvDxi>());
		auto & dcsdx = dcdx->back();
		for (auto &d : dadx)
			dcsdx.push_back((c * n) * d);
	}
	return res;
}
PhiPsiCodes::PhiPsiCodes(const std::vector<NSPgeometry::XYZ> &crd,
		std::vector<BSInChain>::const_iterator bs, bool Nterm, bool Cterm) {
	using namespace NSPgeometry;
	if (!Nterm) {
		std::vector<XYZ> dphi;
		phi = torsion(crd[(bs - 1)->cid], crd[bs->nid],
				crd[bs->caid], crd[bs->cid], &dphi);
		dphidx.push_back(std::make_pair((bs - 1)->cid, dphi[0]));
		dphidx.push_back(std::make_pair(bs->nid, dphi[1]));
		dphidx.push_back(std::make_pair(bs->caid, dphi[2]));
		dphidx.push_back(std::make_pair(bs->cid, dphi[3]));
		phicodes = gettriangcodes(phi, dphidx, &dphicodesdx);
	}
	if (!Cterm) {
		std::vector<XYZ> dpsi;
		psi = torsion(crd[bs->nid], crd[bs->caid],
				crd[bs->cid], crd[(bs + 1)->nid], &dpsi);
		dpsidx.push_back(std::make_pair(bs->nid, dpsi[0]));
		dpsidx.push_back(std::make_pair(bs->caid, dpsi[1]));
		dpsidx.push_back(std::make_pair(bs->cid, dpsi[2]));
		dpsidx.push_back(std::make_pair((bs + 1)->nid, dpsi[3]));
		psicodes = gettriangcodes(psi, dpsidx, &dpsicodesdx);
	}
}


void IntrctMol::calcphipsicodes(bool updateonly) const{
	recollectcrds_all();
	if(!updateonly) {
		phipsicodes.clear();
		for (auto &bss:bsinchains){
			phipsicodes.push_back(makephipsicodes(crds_all_,bss));
		}
	} else {
		for(int cid=0;cid<bsinchains.size();++cid){
#pragma omp parallel for schedule(dynamic,1)
			for(int pid=0;pid<bsinchains[cid].size();++pid){
				if(windowactive(bsinchains[cid],pid,1))
					phipsicodes[cid][pid]=PhiPsiCodes(crds_all_,
							bsinchains[cid].cbegin()+pid,pid == 0,
							pid==bsinchains[cid].size() - 1);
			}
		}
	}
}
void IntrctMol::calcsscodes(bool updateonly) const{
	if(phipsicodes.empty()) calcphipsicodes();
	if(!updateonly) {
		sscodes.clear();
		for (auto &ppc:phipsicodes) {
			sscodes.push_back(estimatess(ppc));
		}
	} else {
		for(int ic=0;ic<bsinchains.size();++ic){
#pragma omp parallel for schedule(dynamic,1)
			for (int i = 0; i < bsinchains[ic].size(); ++i) {
				if (!windowactive(bsinchains[ic],i,3))
					continue;
				auto & sc = sscodes[ic].at(i);
				sc.p3 = NN_SSTerm().probabilities(phipsicodes[ic], i,
						&(sc.dp3dx));
				sc.ssid = NN_SSTerm::sstype(sc.p3);
			}
		}
	}
}
std::vector<std::string> IntrctMol::sscodestrings() const{
	assert(offsetsok_);
	if(sscodes.empty()) calcsscodes();
	std::vector<std::string> res;
	for( auto &c:sscodes){
		res.push_back(std::string(c.size(),' '));
		int idx=0;
		for(auto &code:c){
             res.back()[idx++]=code.charcode();
		}
	}
	return res;
}

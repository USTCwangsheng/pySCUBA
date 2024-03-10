/*
 * intrctblock.cpp
 *
 *  Created on: 2019年12月1日
 *      Author: hyiu
 */
#include "iblock/blcktopo.h"
#include "iblock/intrctblck.h"
#include "iblock/backboneff.h"
#include "iblock/vsctypes.h"
#include "proteinrep/pdbrecord.h"
using namespace NSPintrct;
NSPgeometry::XYZ IntrctBlck::getrepcrd() const {
	if (isaminoacid) {
		return getcrd(topo->atomseqs.at("CA"));
	} else {
		return NSPgeometry::center(crds);
	}
}
NSPproteinrep::BackBoneSite IntrctBlck::getbackbonesite() const {
	static const std::string chainidchar { "ABCDEFGHIJK" };
	NSPproteinrep::BackBoneSite bs;
	std::vector<NSPgeometry::XYZ> bcrd;
	bcrd.push_back(10.0 * getcrd("N"));
	bcrd.push_back(10.0 * getcrd("CA"));
	bcrd.push_back(10. * getcrd("C"));
	bcrd.push_back(10.0 * getcrd("O"));
	bs.changecrd(bcrd);
	bs.resname = topo->blckname;
	bs.chainid = chainidchar[idx2d_.first];
	bs.pdbid = "XXXX";
	bs.resseq = idx2d_.second + 1;
	bs.resid = idx2d_.second + 1;
	bs.data_[NSPproteinrep::BackBoneSite::PHI] = 360.0;
	bs.data_[NSPproteinrep::BackBoneSite::PSI] = 360.0;
	bs.data_[NSPproteinrep::BackBoneSite::OMIGA] = 180.0;
	return bs;
}
void IntrctBlck::changesccrds(const std::vector<NSPgeometry::XYZ> &sccrds) {
	assert(sccrds.size() == topo->natoms - 4);
	std::copy(sccrds.begin(), sccrds.end(), crds.begin() + 2);
}
std::vector<IntrctBlck::AtomNbList> IntrctBlck::getneighbors(
		const IntrctBlck &blck2, const IntrctPara &param,
		std::array<int, 3> &nnb) const {
	std::vector<AtomNbList> res(3, AtomNbList(topo->natoms, std::set<int>()));
	AtomNbList & nblist_M = res[0];
	AtomNbList & nblist_S = res[1];
	AtomNbList &nblist_L = res[2];
	bool proexcl = prepro
			&& (blck2.getidx2d().second == idx2d_.second + 1
					&& blck2.getidx2d().first == idx2d_.first);
	int ibegin = 0, iend = topo->natoms;
	if (activemod == SIDECHAIN && blck2.activemod == INACTIVE) {
		ibegin = 2;
		iend -= 2;
	}
	bool protein = (this->topo->atomips[0].isprotein)
			&& (blck2.topo->atomips[0].isprotein);
	for (int i = ibegin; i < iend; ++i) {
		int jbegin = 0, jend = blck2.topo->natoms;
		if (blck2.activemod == SIDECHAIN && activemod == INACTIVE) {
			jbegin = 2;
			jend -= 2;
		}  //skip checking main chain atoms
		for (int j = jbegin; j < jend; ++j) {
			AtomNbList *l = &nblist_M;
			int *nb = &(nnb[0]);
			if (protein) {
				if (this->topo->atomips[i].issidechain
						|| blck2.topo->atomips[j].issidechain) {
					l = &nblist_S;
					nb = &(nnb[1]);
				}
			} else {
				AtomNbList *l = &nblist_L;
				nb = &(nnb[2]);
			}
			double d2 = crds[i].squaredDistance(blck2.getcrd(j));
			if (d2 <= param.max_cutoff2) {
				if (this->topo->excld[i].find(
						torelative_index(
								AtomIdx3D { blck2.idx2d_.first,
										blck2.idx2d_.second, j }))
						!= this->topo->excld[i].end())
					continue;
				if (proexcl) {
					if (j == topo->CD_PRO) {
						if (i == 1 || i >= topo->natoms - 2)
							continue;
					}
				}
				if (!protein) {
					if (lexcludedatoms.find(j) != lexcludedatoms.end()
							|| blck2.lexcludedatoms.find(i)
									!= blck2.lexcludedatoms.end())
						continue;
				}
				l->at(i).insert(j);
				++(*nb);
			}
		}
	}
	return res;
}
void IntrctBlck::mksteric_nblist(const IntrctPara &param) const {
	nblist_M.clear();
	nblist_S.clear();
	nblist_L.clear();
	if (activemod != INACTIVE) {
//		if(activemod ==ALL) insert_nblist_M(*this,topo->nblist_internal_M);
//		insert_nblist_S(*this,topo->nblist_internal_S);
	}

	for (int i1 = idx2d_.first; i1 < molsystm_->sized1(); ++i1) { //current and later chains
		int i2_begin = (i1 == idx2d_.first) ? idx2d_.second + 1 : 0;
		for (int i2 = i2_begin; i2 < molsystm_->sized2(i1); ++i2) { //subsequent blocks
			if (activemod == INACTIVE
					&& (*molsystm_)(i1, i2).activemod == INACTIVE)
				continue; //skip inactive-inactive blocks;
			if (regioncode * (*molsystm_)(i1, i2).regioncode < 0)
				continue;
			std::array<int, 3> nnb { 0, 0, 0 };
			auto neighbors = getneighbors((*molsystm_)(i1, i2), param, nnb);
			if (nnb[0] > 0)
				insert_nblist_M((*molsystm_)(i1, i2), neighbors[0]);
			if (nnb[1] > 0)
				insert_nblist_S((*molsystm_)(i1, i2), neighbors[1]);
			if (nnb[2] > 0)
				insert_nblist_L((*molsystm_)(i1, i2), neighbors[2]);
		}
	}
}
using namespace NSPproteinrep;
std::vector<PdbRecord> IntrctBlck::topdbrecords() const {
	const std::vector<char> chainids { 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
			'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U',
			'V', 'W', 'X', 'Y', 'Z' };
	std::vector<PdbRecord> records;
	int elementsymbolsize = 1;
	for (int i = 0; i < natoms(); ++i) {
		PdbRecord record;
		record.label = "ATOM";
		record.atomname = atomname(i);
		//get atomname from key, and assign the atomname to record.atomname
		record.namesymbol = record.atomname.substr(0, elementsymbolsize);
		if (elementsymbolsize > 1) {
			record.elementname[0] = record.namesymbol[0];
			record.elementname[1] = record.namesymbol[1];
		} else {
			record.elementname[1] = record.namesymbol[0];
		}
		record.namemodifier = record.atomname.substr(elementsymbolsize);
		record.residuename = resname();
		record.atomid = aoffset + i;
		record.residueid = idx2d_.second + 1;
		record.x = crds[i].x_ * 10.0;
		record.y = crds[i].y_ * 10.0;
		record.z = crds[i].z_ * 10.0;
		int chainid = idx2d_.first;
		assert(chainid < 26);
		record.chainid = chainids[chainid];
		records.push_back(PdbRecord(record.toString()));
	}
	return records;
}


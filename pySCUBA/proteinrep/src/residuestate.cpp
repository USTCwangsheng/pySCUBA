/*
 * residuestate.cpp
 *
 *  Created on: 2020年12月29日
 *      Author: hyiu
 */

#include "proteinrep/residuestate.h"
using namespace NSPproteinrep;
using namespace NSPdesignseq;
std::vector<std::vector<ResidueState>> NSPproteinrep::residuestates(
		NSPdesignseq::PDB &pdb){
	NSPdesignseq::StructureInfo strinfo(&pdb);
	strinfo.updateTorsion();
	strinfo.updateSecondaryStructure();
	SasaPSD sasapsd;
	strinfo.updateSAI(&sasapsd);
	std::vector<ProteinChain *> &pchains=pdb.getChains();
	int nchains=pchains.size();
	std::vector<std::vector<ResidueState>> result;
	for(int c=0;c<nchains;++c){
		int nres=pchains[c]->getResList().size();
		result.push_back(std::vector<ResidueState>(nres));
	}
	int seqid=0;
	for(int c=0;c<nchains;++c){
		for(auto &state:result[c]){
			BackBoneState &bs=state.backbonestate;
			SideChainState &scs=state.sidechainstate;
			bs.phi=strinfo.getPhi(seqid);
			bs.psi=strinfo.getPsi(seqid);
			bs.omiga=strinfo.getOmg(seqid);
			bs.SSState=strinfo.getSS(seqid);
			bs.sai=strinfo.getSai(seqid);
			Residue *residue=strinfo.getResidue(seqid);
			scs.resiudetype=residue->getType();
			scs.torsions=sidechaintorsions(residue);
			seqid++;
		}
	}
	return result;
}

std::vector<double> NSPproteinrep::kai_diffs(const SideChainState &s1, const SideChainState & s2,
		bool loose_degenerate){
	std::vector<double> diffs;
	int nkai=s1.torsions.size()<s2.torsions.size()? s1.torsions.size():s2.torsions.size();
	for(int i=0;i<nkai;++i){
		double diff=s2.torsions.at(i)-s1.torsions.at(i);
		if(diff>180.0) diff-=360.0;
		else if(diff<=-180.0) diff+=360.0;
		if(diff<-90.0 || diff>90.0){
			int (*degen)(const std::string &);
			if (loose_degenerate) degen= &(NSPproteinrep::degeneratedkai_loose);
			else degen=&(NSPproteinrep::degeneratedkai_strict);
			if(i==degen(s1.resiudetype) || i==degen(s2.resiudetype)){
				if(diff<-90.0) diff+=180.0;
				else diff -=180.0;
			}
		}
		diffs.push_back(diff);
	}
	return diffs;
}
int NSPproteinrep::degeneratedkai_loose(const std::string &resname){
	if(resname=="ASP" || resname=="ASN" ||
			resname=="HIS"||
			resname=="PHE" ||
			resname=="TYR") return 1;
	else if(resname=="GLU" ||resname=="GLN") return 2;
	else return -1;
}
int NSPproteinrep::degeneratedkai_strict(const std::string &resname){
	if(resname=="ASP" ||
			resname=="PHE" ||
			resname=="TYR") return 1;
	else if(resname=="GLU") return 2;
	else return -1;
}
#include "iblock/vsctypes.h"
using namespace NSPintrct;
void NSPproteinrep::checkatomnames(NSPdesignseq::PDB &pdb){
	auto residues=pdb.getResList();
	for(auto r:residues){
		const VSCType & vsc=VSCType::getVSCType(r->getType());
		for(auto &anm:vsc.atomnames) {
			if(r->hasAtom(anm)) continue;
			if(anm.size()>2){
				std::string anm1=anm.substr(0,2);
				if(r->hasAtom(anm1)){
					r->renameatom(anm1,anm);
				}
			}
		}
		if((!r->hasAtom("O")) && r->hasAtom("OT1")){
						r->renameatom("OT1","O");
		}
	}
}
std::vector<double> NSPproteinrep::sidechaintorsions(NSPdesignseq::Residue *residue){
	const VSCType & vsc=VSCType::getVSCType(residue->getType());
	const std::vector<int> & rotameratoms=vsc.rotameratoms;
	std::vector<double> result;
	const double rad=180.0/3.14159265;
	for(int l:rotameratoms){
		std::vector<NSPgeometry::XYZ> xyz;
		if(!residue->hasAtom(vsc.atomnames.at(l))) return result;
		xyz.push_back(residue->getAtom(vsc.atomnames.at(l))->getCoord());
		for (int m=0;m<3;++m){
			int a=vsc.internalcrds[l][m].first;
			std::string aname;
			if(a==0) aname="N";
			else if(a==1) aname="CA";
			else aname=vsc.atomnames[a-2];
			xyz.push_back(residue->getAtom(aname)->getCoord());
		}
		double t=NSPgeometry::torsion(xyz[0],xyz[1],xyz[2],xyz[3])*rad;
		result.push_back(t);
	}
	return result;
}

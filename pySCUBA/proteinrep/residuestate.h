/*
 * residuestate.h
 *
 *  Created on: 2020年12月29日
 *      Author: hyiu
 */

#ifndef PROTEINREP_RESIDUESTATE_H_
#define PROTEINREP_RESIDUESTATE_H_
#include "designseq/StructureInfo.h"
namespace NSPproteinrep{
struct BackBoneState{
	double phi{360.0};
	double psi{360.0};
	double omiga{180.0};
	double sai{-1.0};
	int phipsistate{-1};
	char SSState{' '};
	int SASAState{-1};
};
struct SideChainState{
	std::string resiudetype{"UNK"};
	int rotamerstate{-1};
	std::vector<double> torsions;
};
struct ResidueState{
	BackBoneState backbonestate;
	SideChainState sidechainstate;
};
void checkatomnames(NSPdesignseq::PDB &pdb);
std::vector<std::vector<ResidueState>> residuestates(NSPdesignseq::PDB &pdb);
inline std::vector<std::vector<ResidueState>> residuestates(std::istream &is){
	NSPdesignseq::PDB pdb(is,"NNNN");
	checkatomnames(pdb);
	return residuestates(pdb);
}
std::vector<double> sidechaintorsions(NSPdesignseq::Residue *residue);
int degeneratedkai_loose(const std::string &resname);
int degeneratedkai_strict(const std::string &resname);
std::vector<double> kai_diffs(const SideChainState &s1,const SideChainState &s2,
		bool loose_degenerate=true);
}

#endif /* PROTEINREP_RESIDUESTATE_H_ */

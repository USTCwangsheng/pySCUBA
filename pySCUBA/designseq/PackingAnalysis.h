/*
 * PackingAnalysis.h
 *
 *  Created on: 2018Äê4ÔÂ24ÈÕ
 *      Author: notxp
 */

#ifndef DESIGNSEQ_PACKINGANALYSIS_H_
#define DESIGNSEQ_PACKINGANALYSIS_H_


#include <string>
#include <vector>
#include "designseq/ProteinRep.h"
#include "designseq/StructureInfo.h"


namespace NSPdesignseq{


using namespace std;
using namespace NSPdesignseq;
using namespace NSPgeometry;

int atomicRMS( PDB* pdbA,  PDB* pdbB,  AtomLib* atLib, vector<float>& distanceList, vector<float>& saiList, vector<string>& typeList);

int atomicRMS( PDB* pdbA,  PDB* pdbB,  AtomLib* atLib, vector<float>& distanceList, string triName);

void chi1Diff( PDB* pdbA,  PDB* pdbB,AtomLib* atLib, vector<float>& chi1DiffList, vector<float>& saiList, vector<string>& typeList);

void chi1Diff( PDB* pdbA,  PDB* pdbB,AtomLib* atLib,  vector<float>& chi1DiffList, string triName);

void chi2Diff( PDB* pdbA,  PDB* pdbB, AtomLib* atLib,  vector<float>& chi2DiffList, vector<float>& saiList, vector<string>& typeList);

void disulfideBondDiff( PDB* pdbA,  PDB* pdbB,  vector<string>& inBoth, vector<string>& onlyA, vector<string>& onlyB);

//void hydrogenBondDiff( PDB* pdbA,  PDB* pdbB,  vector<string>& inBoth, vector<string>& onlyA, vector<string>& onlyB);

//void hydrogenBondCount(PDB* pdb, AtomLib* atLib,  int* counts, set<string>& hbs);

float squareAverage(vector<float>& list);

}


#endif /* ANALYSIS_SIDECHAINCOMPARISION_H_ */

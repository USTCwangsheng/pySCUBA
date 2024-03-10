/*
 * ScoringTemplate.h
 *
 *  Created on: Dec 20, 2018
 *      Author: s2982206
 */

#ifndef DESIGNSEQ_SCORINGTEMPLATE_H_
#define DESIGNSEQ_SCORINGTEMPLATE_H_
#include <time.h>
#include "designseq/ProteinRep.h"
#include "designseq/Conformer.h"
#include "designseq/S1EnergyTable.h"
#include "designseq/S2MatrixFinder.h"
#include "designseq/SeqProfile.h"
#include "designseq/Rotamer.h"
#include "designseq/RotamerLib.h"
#include "designseq/EnergyArray.h"
#include "designseq/EnergyMatrix.h"
#include "designseq/DesignParameters.h"
#include "designseq/StructureInfo.h"
#include "designseq/AtomicEnergyCalcular.h"
#include "designseq/AtomLib.h"
#include "designseq/AADesignTemplate.h"
#include "backbone/backbonesite.h"

namespace NSPdesignseq {
using namespace NSPproteinrep;

class ScoringTemplate {


public:
    vector<EnergyArray*> s1EAList;
    vector<EnergyArray*> refEAList;
    vector<EnergyArray*> rotEAList;
    vector<EnergyMatrix*> s2EMList;
    vector<EnergyMatrix*> packEMList;

    vector<Residue*> resList;
    vector<BackboneSitesPair*> bsPairs;
    vector<vector<int>*> involvedPairs;
    vector<RotamerGroup*> rotGroups;
    vector<BackBoneSite*> bsList;
    int resNum;
    map<string, int> chainIDResIDToIndex;
    vector<vector<int>*> neighborResidues;

    vector<AAProbabilityArray*> paList;
    vector<AAScoreMatrix*> smList;

    RotamerLib* rotLibA;
    RotamerLib* rotLibB;

    AtomLib atLib;
    AtomicEnergyCalcular* aec;
    ResName rn;
    DesignParameters* dp;

    S1EnergyTable* s1ET;



    ScoringTemplate(PDB* targetPDB, DesignParameters* dp, S1EnergyTable& s1Etable, S2MatrixFinder& s2Etable);
    void loadSingleResidueEnergy();
    void loadPairwiseEnergy();
    void updateResPairs();
    void loadS1S2(S1EnergyTable& s1Etable, S2MatrixFinder& s2Etable);
    bool contact(int posA, int posB, float cutoff);
    void setRotamerGroups();
    float clashScore(vector<Atom*>& listA, vector<Atom*>& listB, string typeB);
    void deleteSc_BB_ClashedRotamers(float cutoff);

    virtual ~ScoringTemplate();
};

float pairEnergyScoring(Conformer* confA, Conformer* confB,int seqSep, float wt, AtomicEnergyCalcular* aec);


} /* namespace NSPtest */

#endif /* DESIGNSEQ_SCORINGTEMPLATE_H_ */

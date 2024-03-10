/*
 * PackingTemplate.h
 *
 *  Created on: 2018Äê4ÔÂ16ÈÕ
 *      Author: notxp
 */

#ifndef DESIGNSEQ_PACKINGTEMPLATE_H_
#define DESIGNSEQ_PACKINGTEMPLATE_H_
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
#include "designseq/DesignTemplate.h"


namespace NSPdesignseq {
using namespace std;

class PackingTemplate : public EnergyCalculatorTemplate{
private:


    PDB* pdb;
    vector<Residue*> resList;

//    vector<BackBoneSite*> bsList;
//    vector<BackboneSitesPair*> bsPairs;

//    vector<vector<int>*> involvedPairs;
    vector<vector<int>*> neighborResidues;


//    vector<EnergyArray*> eaList;
//    vector<EnergyMatrix*> emList;

//    vector<RotamerGroup*> rotGroups;

    PhipsiLib* ppLib;

    RotamerLib* rotLibA;
    RotamerLib* rotLibB;

    AtomLib* atLib;
    ResName rn;


    DesignParameters* para;

public:
    AtomicEnergyCalcular* ec;
    int resNum;
    PackingTemplate(PDB* pdb, DesignParameters* para);


    void updateRotamerGroups();
    void updateResPairs();
    bool contact(int posA, int posB, float cutoff);

    void updateEAEM();

    void printPositionChoiceNum();

    void printDetail();
    void testUbq();

    virtual ~PackingTemplate();
};


float pairEnergyPack(Conformer* confA, Conformer* confB,int seqSep, float wt, AtomicEnergyCalcular* aec);
float pairHBEnergyPack(Conformer* confA, Conformer* confB, int index, AtomicEnergyCalcular* aec);

} /* namespace packing */

#endif /* DESIGNSEQ_PACKINGTEMPLATE_H_ */

/*
 * SingleSitePackingTemplate.h
 *
 *  Created on: Sep 21, 2018
 *      Author: s2982206
 */

#ifndef DESIGNSEQ_SINGLESITEPACKINGTEMPLATE_H_
#define DESIGNSEQ_SINGLESITEPACKINGTEMPLATE_H_
#include <time.h>
#include "designseq/DesignTemplate.h"
#include "designseq/PackingTemplate.h"
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

class SingleSitePackingTemplate: public EnergyCalculatorTemplate {
public:
    int resNum;
    map<string,int> chainIDResIDToIndex;
    vector<vector<int>*> neighborResidues;
    set<int> fixedPos;

    vector<Residue*> resList;
    vector<Rotamer*> natRotList;

    RotamerLib* rotLibA;
    RotamerLib* rotLibB;

    AtomLib atLib;
    AtomicEnergyCalcular aec;
    ResName rn;
    DesignParameters* dp;
    int mutPosition;
    string mutType;

    SingleSitePackingTemplate(ProteinChain* pc, DesignParameters* dp, char chainID, string resID, string triName);

    void nativeAARotamerGroup();
    void nonContactRotamer();
    void loadSingleResidueEnergy();
    void loadPairwiseEnergy();
    void updateResPairs();
    void printChoice();

    bool contact(int posA, int posB, float cutoff);

    virtual ~SingleSitePackingTemplate();
};

}

#endif /* DESIGNSEQ_SINGLESITEPACKINGTEMPLATE_H_ */

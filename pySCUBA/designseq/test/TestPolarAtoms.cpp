/*
 * TestPolarAtoms.cpp
 *
 *  Created on: 2018Äê4ÔÂ27ÈÕ
 *      Author: notxp
 */


#include <string>
#include <vector>
#include <iostream>

#include "designseq/Rotamer.h"
#include "designseq/ResName.h"
#include "designseq/RotamerGroup.h"
#include "geometry/xyz.h"
#include "geometry/localframe.h"
#include "designseq/AtomLib.h"
#include "designseq/ProteinRep.h"
#include "designseq/RotamerLib.h"
#include "designseq/Conformer.h"
#include "backbone/backbonesite.h"
#include "designseq/StructureInfo.h"
#include "designseq/DesignParameters.h"
#include "designseq/AtomicEnergyCalcular.h"
#include "designseq/PackingTemplate.h"


using namespace std;
using namespace NSPdesignseq;
using namespace NSPgeometry;


int main(int argc, char** args){


    ResName rn;
    string pdbFile = string(args[1]);

    PDB pdb(pdbFile,"unk");

    ProteinChain* pc = pdb.getFirstChain();
    string natSeq = pc->getSequence();

    vector<Residue*> resList = pc->getResList();


    vector<BackBoneSite*> bsList;
    proteinchain2BackboneSiteList(pc, bsList);
    AtomLib atLib;




    vector<Conformer*> cfList;

    for(int i=0;i<1;i++){
        string triname = bsList[i]->resname;
        Rotamer* rot = resList[i]->natRotamer(&atLib);
        Conformer* cf = new Conformer(rot, bsList[i], &atLib);
        cfList.push_back(cf);
    }

    for(int i=1;i<bsList.size();i++){
        string triname = bsList[i]->resname;
        Rotamer* rot = resList[i]->natRotamer(&atLib);
        Conformer* cf = new Conformer(rot, bsList[i], bsList[i-1], &atLib);
        cfList.push_back(cf);
    }

    DesignParameters dp(2);
    AtomicEnergyCalcular ec(&dp);


    Conformer* cfA = cfList[54];
    Conformer* cfB = cfList[56];
    PolarAtom* pa = cfA->scPolarList[0];
    PolarAtom* pb = cfB->scPolarList[0];
    float e = ec.hbEnergy(pa,pb);
    printf("e= %7.3f\n",e);


}

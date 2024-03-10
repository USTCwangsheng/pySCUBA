/*
 * DSBondEnergy.cpp
 *
 *  Created on: 2018Äê4ÔÂ24ÈÕ
 *      Author: notxp
 */





#include <iostream>
#include <vector>
#include <string>
#include "designseq/ProteinRep.h"
#include "designseq/StructureInfo.h"
#include "designseq/S1EnergyTable.h"
#include "geometry/xyz.h"
#include "designseq/AtomLib.h"
#include "designseq/StructureInfo.h"
#include "designseq/S2MatrixFinder.h"
#include "designseq/AtomicEnergyCalcular.h"
#include "designseq/DesignTemplate.h"
#include "designseq/AAProbabilityArray.h"

using namespace std;
using namespace NSPdesignseq;
using namespace NSPgeometry;

int main(int argc, char** args){

    /*
     * input is PDB structure
     */
    if(argc != 2)
    {
        cout << "testScoring $PDBFile" << endl;
    }
    string s(args[1]);
    string pdbID = "unk";
    PDB pdb(s, pdbID);
    ProteinChain* pc = pdb.getFirstChain();

    /*
     * calculate backbone torsion angle, secondary structure, and SASA
     */
    StructureInfo strInfo(pc);
    int resNum = pc->getChainLength();
    SasaPSD sasaPSD;
    strInfo.updateTorsion();
    strInfo.updateSecondaryStructure();
    strInfo.updateSAI(&sasaPSD);





    vector<BackBoneSite*> bsList;
    proteinchain2BackboneSiteList(pc, bsList);


    /*
     * calculate atomic score
     */

    double totVdwEnergy = 0;

    AtomicEnergyCalcular ec;
    DesignParameters dp;
    AtomLib atLib;
    vector<Rotamer*> rotList;
    for(int i=0;i<resNum;i++){
        Residue* res = pc->getResList().at(i);
        rotList.push_back(res->natRotamer(&atLib));
    }

    for(int i=0;i<resNum;i++){
        string resIDA = pc->getResList().at(i)->resID;
        BackBoneSite* bsA = (bsList.at(i));
        Rotamer* rotA = rotList.at(i);
        string triNameA = rotA->triName;
        if(triNameA != "CYS") continue;

        Conformer* confA = new Conformer(rotA, bsA, &atLib);
        double saiA = strInfo.getSai(i);

        for(int j=i+2;j<resNum;j++){
            string resIDB = pc->getResList().at(j)->resID;

            double saiB = strInfo.getSai(j);
            double vdwWT = dp.getVdwWeight(0.5*saiA + 0.5*saiB);
            BackBoneSite* bsB = (bsList.at(j));
            Rotamer* rotB = rotList.at(j);
            string triNameB = rotB->triName;
            if(triNameB != "CYS") continue;
            Conformer* confB = new Conformer(rotB, bsB, &atLib);

            XYZ* cbI = confA->scCoordList[0];
            XYZ* sgI = confA->scCoordList[1];
            XYZ* cbJ = confB->scCoordList[0];
            XYZ* sgJ = confB->scCoordList[1];
            float dSS = sgI->distance(*sgJ);
            float ang1 = NSPdesignseq::angle(*cbI, *sgI, *sgJ);
            float ang2 = NSPdesignseq::angle(*sgI, *sgJ, *cbJ);
            float dihed = dihedral(*cbI, *sgI, *sgJ, *cbJ);

            float eSS = (dSS-2.055)*(dSS-2.055);
            eSS += 0.00063*(ang1-105.2)*(ang1-105.2);
            eSS += 0.00063*(ang2-105.2)*(ang2-105.2);
            eSS += 0.00018*(abs(dihed)-89.7)*(abs(dihed)-89.7);
            if(eSS < 5.0)
            {
                printf("Energy: %7.4f\n",eSS);
            }
        }
    }





}

/*
 * TestSingleMutation.cpp
 *
 *  Created on: Sep 21, 2018
 *      Author: s2982206
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
#include "designseq/SeqMinimize.h"
#include "designseq/AAProbabilityArray.h"
#include "designseq/SingleSitePackingTemplate.h"

using namespace std;
using namespace NSPdesignseq;
using namespace NSPgeometry;



int main(int argc, char** args){

    /*
     * input is PDB structure
     */

    clock_t start = clock();

    if(argc != 6)
    {
        cout << "singleMut $PDBFile $ChainID $ResID $MutType" << endl;
        exit(0);
    }
    string s(args[1]);
    string pdbID = "unkn";

    if(s.length() > 7) {
        int x = 0;
        int y = 0;
        for(int i=0;i<s.length();i++) {
            if(s[i] == '/')
                x = i+1;
            if(s[i] == '.')
                y = i;
        }
        pdbID = s.substr(x,y-x);
    }


    PDB pdb(s, pdbID);

    DesignParameters dp;

    char chainID = args[2][0];
    string resid(args[3]);
    string mutType(args[4]);
    string output(args[5]);

    ProteinChain* pc = pdb.getChain(chainID);
    Residue* res = pc->getResidue(resid);
    string natType = res->triName;

    SingleSitePackingTemplate* st1 = new SingleSitePackingTemplate(pc, &dp, chainID, resid, mutType);
    int seqLen = st1->resNum;
    RotSequence unit1(seqLen);
    DesignMC mc1(st1);
    mc1.packing2(&unit1);
    double packEnergy1 = unit1.totEnergy(st1);
    cout << "packEnergy: " << packEnergy1 << endl;


    SingleSitePackingTemplate* st2 = new SingleSitePackingTemplate(pc, &dp, chainID, resid, natType);

    RotSequence unit2(seqLen);
    DesignMC mc2(st2);
    mc2.packing2(&unit2);
    double packEnergy2 = unit2.totEnergy(st2);

    double packEnergy = packEnergy1 - packEnergy2;

    vector<Residue*> tmpResList = pc->getResList();
    vector<Residue*> resList;
    int index = 0;

    for(unsigned int i=0;i<tmpResList.size();i++)
    {
        Residue* res = tmpResList.at(i);
        if(res->hasThreeCoreAtoms()&& res->hasAtom("O"))
        {
            resList.push_back(res);
        }
    }

    StructureInfo si(pc);
    si.updateBBHBonds();

    vector<BackBoneSite*> bsList;
    proteinchain2BackboneSiteList(pc, bsList);

    int resNum = bsList.size();


    if(res == NULL){
        cerr << "invalid resid: " << resid << endl;
        exit(0);
    }

    int resSeqID = -1 ;
    for(int i=0;i<resNum;i++){
        if(resList[i]->resID == resid)
            resSeqID = i;
    }
    if(resSeqID == -1){
        cerr << "residue backbone information insufficient!" << endl;
        abort();
    }



    ResName rn;
    /*
     * calculate s1 score
     */

    /*
    S1EnergyTable s1ET;

    AAProbabilityArray pa;
    s1ET.getS1(*bsList[resSeqID], &pa);
    int intNameNat = rn.triToInt(res->triName);
    int intNameMut = rn.triToInt(mutType);
    if(intNameMut > 19){
        cerr << "incorrect residue type, please input the mutant type in three letter code" << endl;
        exit(0);
    }

    float s1Nat = pa.getScore(intNameNat);
    float s1Mut = pa.getScore(intNameMut);
    */

//    printf("s1 nat: %7.3f mut: %7.3f\n",s1Nat, s1Mut);

    /*
     * calculate S2 score
     */

    /*
    S2MatrixFinder s2ET;
    float s2Nat = 0;
    float s2Mut = 0;

    BackBoneSite* bsTarget = bsList.at(resSeqID);
    for(int i=0;i<resSeqID;i++){
        BackBoneSite* bsA = bsList[i];
        int intNameA = rn.triToInt(bsA->resname);
        float cbDistance = bsA->cbcrd().distance(bsTarget->cbcrd());
        if(cbDistance > 8.0) continue;
        BackboneSitesPair pair(bsA, bsTarget);
        AAScoreMatrix sm;
        s2ET.getSM(&pair,&sm);
        s2Nat += sm.getValue(intNameA, intNameNat);
        s2Mut += sm.getValue(intNameA, intNameMut);

//        printf("s2 pair: %3d %3d nat: %7.3f mut: %7.3f\n",bsA->resid,bsTarget->resid,sm.getValue(intNameA, intNameNat),  sm.getValue(intNameA, intNameMut) );
    }

    for(int i=resSeqID+1;i<resNum;i++){
        BackBoneSite* bsB = bsList[i];
        int intNameB = rn.triToInt(bsB->resname);
        float cbDistance = bsB->cbcrd().distance(bsTarget->cbcrd());
        if(cbDistance > 8.0) continue;
        BackboneSitesPair pair(bsTarget, bsB);
        AAScoreMatrix sm;
        s2ET.getSM(&pair,&sm);
        s2Nat += sm.getValue(intNameNat, intNameB);
        s2Mut += sm.getValue(intNameNat, intNameB);
//        printf("s2 pair: %3d %3d nat: %7.3f mut: %7.3f\n",bsTarget->resid,bsB->resid, sm.getValue(intNameNat, intNameB),  sm.getValue(intNameNat, intNameB));

    }

    double saiA = bsTarget->data_[3];
    printf("%4s %c %4s %3s->%3s SS: %c SAI: %5.3f S1: %6.3f S2: %6.3f PACK: %6.3f\n", pdbID.c_str(), chainID, resid.c_str(), res->triName.c_str(), mutType.c_str(), bsTarget->sscode, saiA, s1Mut-s1Nat, s2Mut-s2Nat, packEnergy);
    */



    clock_t end = clock();
    cout << "time: " << (float)(end-start)/CLOCKS_PER_SEC << endl;
}



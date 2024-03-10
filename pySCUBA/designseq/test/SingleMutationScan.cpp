/*
 * SingleMutationScan.cpp
 *
 *  Created on: 2018Äê4ÔÂ22ÈÕ
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
#include "designseq/CmdArgs.h"
#include "designseq/SeqMinimize.h"
#include "designseq/SingleSitePackingTemplate.h"

using namespace std;
using namespace NSPdesignseq;
using namespace NSPgeometry;


vector<double> natPackEnergyAndHbondEnergy(ProteinChain* pc, string resID, vector<BackBoneSite*> bsList, int resSeqID){
    BackBoneSite* bsTarget = bsList[resSeqID];
    DesignParameters dp;

    SingleSitePackingTemplate* stNat = new SingleSitePackingTemplate(pc, &dp, pc->getChainID(), resID, bsTarget->resname);
    RotSequence* unitNat = new RotSequence(stNat->resNum);
    DesignMC* mcNat = new DesignMC(stNat);
    mcNat->packing3(unitNat);
    double energyNatPack = unitNat->totEnergy(stNat);
    double natHB = unitNat->resInvolvedHbondEnergy(resSeqID, stNat);
    vector<double> list;
    list.push_back(energyNatPack);
    list.push_back(natHB);
    delete stNat;
    delete unitNat;
    delete mcNat;
    return list;
}

string getMutationInfo(ProteinChain* pc, string resID, vector<BackBoneSite*> bsList, int resSeqID, string mutType, S1EnergyTable* s1ET, S2MatrixFinder* s2ET) {

    ResName rn;
    BackBoneSite* bsTarget = bsList[resSeqID];
    int resNum = bsList.size();

    AAProbabilityArray pa;
    s1ET->getS1(*bsList[resSeqID], &pa);
    int intNameNat = rn.triToInt(bsTarget->resname);
    int intNameMut = rn.triToInt(mutType);
    if(intNameMut > 19){
        cerr << "invalid residue type, please input the mutant type in three letter code" << endl;
        exit(0);
    }


    float s1Nat = pa.getScore(intNameNat);
    float s1Mut = pa.getScore(intNameMut);
    float s2Nat = 0;
    float s2Mut = 0;


    for(int i=0;i<resSeqID;i++){
        BackBoneSite* bsA = bsList[i];
        int intNameA = rn.triToInt(bsA->resname);
        float cbDistance = bsA->cbcrd().distance(bsTarget->cbcrd());
        if(cbDistance > 8.0) continue;
        BackboneSitesPair pair(bsA, bsTarget);

        AAScoreMatrix sm;
        s2ET->getSM(&pair,&sm);
        s2Nat += sm.getValue(intNameA, intNameNat);
        s2Mut += sm.getValue(intNameA, intNameMut);
    }

    for(int i=resSeqID+1;i<resNum;i++){
        BackBoneSite* bsB = bsList[i];
        int intNameB = rn.triToInt(bsB->resname);
        float cbDistance = bsB->cbcrd().distance(bsTarget->cbcrd());
        if(cbDistance > 8.0) continue;
        BackboneSitesPair pair(bsTarget, bsB);

        AAScoreMatrix sm;
        s2ET->getSM(&pair,&sm);
        s2Nat += sm.getValue(intNameNat, intNameB);
        s2Mut += sm.getValue(intNameNat, intNameB);

    }

    DesignParameters dp;

    vector<double> natEnergyList = natPackEnergyAndHbondEnergy(pc, resID, bsList, resSeqID);

    double energyNatPack = natEnergyList[0];
    double natHB = natEnergyList[1];

    SingleSitePackingTemplate* stMut = new SingleSitePackingTemplate(pc, &dp, pc->getChainID(), resID, mutType);
    RotSequence* unitMut = new RotSequence(stMut->resNum);
    DesignMC* mcMut = new DesignMC(stMut);
    mcMut->packing3(unitMut);
    double energyMutPack = unitMut->totEnergy(stMut);
    double mutHB = unitMut->resInvolvedHbondEnergy(resSeqID, stMut);

    string pdbID = pc->getPDBID();
    char chainID = pc->getChainID();
    double ePack = energyMutPack - energyNatPack;
    double eHB = mutHB - natHB;
    double saiA = bsTarget->data_[3];
    int chainLen = pc->getChainLength();
    double refNat = 0;
    double refMut = 0;
    if(bsTarget->sscode == 'H') {
        refNat = dp.refH[intNameNat];
        refMut = dp.refH[intNameMut];
    }
    else if(bsTarget->sscode == 'E') {
        refNat = dp.refE[intNameNat];
        refMut = dp.refE[intNameMut];
    }
    else if(bsTarget->sscode == 'C') {
        refNat = dp.refC[intNameNat];
        refMut = dp.refC[intNameMut];
    }

    int helixNum = 0;
    int betaNum = 0;
    int coreNum = 0;
    int surfNum = 0;
    for(int i=0;i<bsList.size();i++){
        double sai = bsList[i]->data_[3];
        char ss = bsList[i]->sscode;
        if(ss == 'H') helixNum += 1;
        if(ss == 'E') betaNum += 1;
        if(sai < 0.33) coreNum += 1;
        if(sai > 0.66) surfNum += 1;
    }

    float pH = 1.0*helixNum/bsList.size();
    float pE = 1.0*betaNum/bsList.size();
    float pCore = 1.0*coreNum/bsList.size();
    float pSurf = 1.0*surfNum/bsList.size();

    char s[200];
    sprintf(s, "%c %4s %3s->%3s SAI: %5.3f S1: %6.3f S2: %6.3f PACK: %7.3f HB: %7.3f", chainID, bsTarget->strresid.c_str(), bsTarget->resname.c_str(), mutType.c_str(), saiA, s1Mut-s1Nat, s2Mut-s2Nat, ePack, eHB);

    delete stMut;
    delete unitMut;
    delete mcMut;

    return string(s);
}


int main(int argc, char** args){

    /*
     * input is PDB structure
     */


    CmdArgs cmd(argc, args);

    if(cmd.specifiedOption("-h"))
    {
        cout << "singleMutScan $PDBFile $ChainID $ResID" << endl;
        cout << "singleMutScan $PDBFile $ChainID" << endl;
        exit(0);
    }


    string s(args[1]);
    string pdbID = "unkn";

    if(s.length() > 7)
        pdbID = s.substr(s.length()-8,4);

    PDB pdb(s, pdbID);

    char chainID = args[2][0];
    ProteinChain* pc = pdb.getChain(chainID);
    if(pc == NULL){
        cerr << "invalid chainID: " << chainID << endl;
        exit(0);
    }

    S1EnergyTable s1ET;
    S2MatrixFinder s2ET;

    StructureInfo si(pc);
    si.updateBBHBonds();
    vector<BackBoneSite*> bsList;
    proteinchain2BackboneSiteList(pc, bsList);
    int resNum = bsList.size();
    string mode = "allSiteScan";
    string targetResID;
    int targetResSeqID = -1;

    if(argc == 3) {

    }
    else if(argc == 4) {
        targetResID = string(args[3]);
        mode = "singleSiteScan";
    }
    else {
        cout << "Usage: " << endl;
        cout << "mutationScan $PDBFile $ChainID $ResID" << endl;
        cout << "mutationScan $PDBFile $ChainID" << endl;
        exit(0);
    }
    ResName rn;

    vector<Residue*> tmpResList = pc->getResList();
    vector<Residue*> resList;
    int index = 0;
    //cout << "tmpList " << tmpResList.size() << endl;
    for(unsigned int i=0;i<tmpResList.size();i++)
    {
        Residue* res = tmpResList.at(i);
        if(res->hasThreeCoreAtoms()&& res->hasAtom("O"))
        {
            resList.push_back(res);
        }
    }

    if(mode == "singleSiteScan") {
        for(unsigned int i=0;i<resList.size();i++) {
            if(resList[i]->resID == targetResID)
                targetResSeqID = i;
        }
        if(targetResSeqID == -1) {
            cout << "invalid resid: " << targetResID << endl;
            exit(0);
        }
        for(int k=0;k<20;k++)
        {
            string mutType = rn.intToTri(k);

            string mutInfo = getMutationInfo(pc, targetResID, bsList, targetResSeqID, mutType, &s1ET, &s2ET);
            cout << mutInfo << endl;
        }
    }
    else {
        for(int i=0;i<resList.size();i++) {
            targetResID = resList[i]->resID;
            targetResSeqID = i;
            for(int k=0;k<20;k++) {
                string mutType = rn.intToTri(k);
                string mutInfo = getMutationInfo(pc, targetResID, bsList, targetResSeqID, mutType, &s1ET, &s2ET);
                cout << mutInfo << endl;
            }
        }
    }



}


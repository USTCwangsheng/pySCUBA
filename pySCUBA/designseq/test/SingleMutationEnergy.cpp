/*
 * SingleMutationEnergy.cpp
 *
 *  Created on: 2018Äê2ÔÂ26ÈÕ
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


    if(argc != 5)
    {
        cout << "singleMutation $PDBFile $ChainID $ResID $MutType" << endl;
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

    char chainID = args[2][0];
    string resid(args[3]);
    string mutType(args[4]);

    ResName rn;

    ProteinChain* pc = pdb.getChain(chainID);
    if(pc == NULL){
        cerr << "invalid chainID: " << chainID << endl;
        exit(0);
    }


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

    Residue* res = pc->getResidue(resid);
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

    int backboneFlexibleIndex = si.getFlexibleIndex(resSeqID);

    /*
     * calculate s1 score
     */
    S1EnergyTable s1ET;

    AAProbabilityArray pa;
    s1ET.getS1(*bsList[resSeqID], &pa);
    int intNameNat = rn.triToInt(res->triName);
    int intNameMut = rn.triToInt(mutType);
    if(intNameMut > 19){
        cerr << "invalid residue type, please input the mutant type in three letter code" << endl;
        exit(0);
    }

    float s1Nat = pa.getRelScore(intNameNat);
    float s1Mut = pa.getRelScore(intNameMut);

    /*
     * calculate S2 score
     */

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
        int sep = pair.seqSep;

        AAScoreMatrix sm;
        s2ET.getSM(&pair,&sm);
        s2Nat += sm.getValue(intNameA, intNameNat);
        s2Mut += sm.getValue(intNameA, intNameMut);
    }

    for(int i=resSeqID+1;i<resNum;i++){
        BackBoneSite* bsB = bsList[i];
        int intNameB = rn.triToInt(bsB->resname);
        float cbDistance = bsB->cbcrd().distance(bsTarget->cbcrd());
        if(cbDistance > 8.0) continue;
        BackboneSitesPair pair(bsTarget, bsB);
        int sep = pair.seqSep;

        AAScoreMatrix sm;
        s2ET.getSM(&pair,&sm);
        s2Nat += sm.getValue(intNameNat, intNameB);
        s2Mut += sm.getValue(intNameNat, intNameB);

    }

    DesignParameters dp;
    SingleSitePackingTemplate* stNat = new SingleSitePackingTemplate(pc, &dp, chainID, resid, res->triName);
    RotSequence* unitNat = new RotSequence(stNat->resNum);
    DesignMC* mcNat = new DesignMC(stNat);
    mcNat->packing3(unitNat);


    double energyNatMin = unitNat->totEnergy(stNat);
    double natHB = unitNat->resInvolvedHbondEnergy(resSeqID, stNat);
    SingleSitePackingTemplate* stMut = new SingleSitePackingTemplate(pc, &dp, chainID, resid, mutType);
    RotSequence* unitMut = new RotSequence(stMut->resNum);
    DesignMC* mcMut = new DesignMC(stMut);
    mcMut->packing3(unitMut);
    double energyMutMin = unitMut->totEnergy(stMut);
    double mutHB = unitMut->resInvolvedHbondEnergy(resSeqID, stMut);

    double eHB = mutHB - natHB;
    double ePack = energyMutMin-energyNatMin;
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
    printf("%c %4s %3s->%3s SAI: %5.3f S1: %6.3f S2: %6.3f PACK: %7.3f HB: %7.3f\n", chainID, resid.c_str(), res->triName.c_str(), mutType.c_str(), saiA, s1Mut-s1Nat, s2Mut-s2Nat, ePack, eHB);

    //printf("%c %4s %3s->%3s SS: %c SAI: %5.3f S1: %6.3f S2: %6.3f PACK: %7.3f HB: %7.3f REF: %5.2f LEN: %3d HELIX: %5.3f SHEET: %5.3f CORE: %5.3f SURF: %5.3f\n", chainID, resid.c_str(), res->triName.c_str(), mutType.c_str(), bsTarget->sscode, saiA, s1Mut-s1Nat, s2Mut-s2Nat, ePack, eHB, refMut-refNat, chainLen, pH, pE, pCore, pSurf);

}


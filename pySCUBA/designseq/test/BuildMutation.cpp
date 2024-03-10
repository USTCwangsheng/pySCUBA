/*
 * BuildMutation.cpp
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
#include "designseq/SeqMinimize.h"
#include "designseq/SingleSitePackingTemplate.h"

using namespace std;
using namespace NSPdesignseq;
using namespace NSPgeometry;


int main(int argc, char** args){

    /*
     * input is PDB structure
     */
    if(argc != 6)
    {
        cout << "buildMutation $PDBFile $ChainID $ResID $MutType $outputFile" << endl;
        exit(0);
    }
    string s(args[1]);
    string pdbID = "unkn";
    string output = args[5];

    if(s.length() > 7)
        pdbID = s.substr(s.length()-8,4);

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
        cerr << "invalid residue: " << resid << endl;
        exit(0);
    }

    int resSeqID = -1 ;
    for(int i=0;i<resNum;i++){
        if(resList[i]->resID == resid)
            resSeqID = i;
    }
    if(resSeqID == -1){
        cerr << "invalid residue: " << chainID << " " << resid << endl;
        abort();
    }


    DesignParameters dp(2);
    SingleSitePackingTemplate* stMut = new SingleSitePackingTemplate(pc, &dp, chainID, resid, mutType);
    RotSequence* unitMut = new RotSequence(stMut->resNum);
    DesignMC* mcMut = new DesignMC(stMut);
    mcMut->packing3(unitMut);
    mcMut->printPDB(unitMut, output);
    delete unitMut;
    delete mcMut;
    delete stMut;

}



/*
 * Test40PDB.cpp
 *
 *  Created on: 2018��1��16��
 *      Author: notxp
 */

#include <iostream>
#include <vector>
#include <string>
#include <time.h>
#include "stdio.h"
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

using namespace std;
using namespace NSPdesignseq;
using namespace NSPgeometry;

int main(int argc, char** args){

    string tag = args[1];


    string parameterFile = "/home/xiongpeng/cpp/designTest/para/"+tag;
    DesignParameters dp(parameterFile);
    S1EnergyTable s1Etable;
    S2MatrixFinder s2Etable;


    string listFile = "/home/xiongpeng/sefTest/list";
    ifstream input;
    input.open(listFile.c_str(),ios::in);
    if (! input.is_open())
    {
        cout << "fail to open file " << listFile << endl;
        exit (1);
    }
    string s;
    vector<string> pdbList;
    while(getline(input,s)){
        pdbList.push_back(s);
    }


    for(int i=0;i<pdbList.size();i++){
        cout << pdbList[i] << endl;
        string pdbFile = "/home/xiongpeng/sefTest/pdb/"+pdbList.at(i)+".pdb";
        string outputFile = "/home/xiongpeng/cpp/designTest/des/" + pdbList[i] + "-" + tag + ".des";
        ofstream output(outputFile, ios::out);
        if(!output.is_open())
        {
            cout << "fail to open file " << outputFile << endl;
            exit(1);
        }

        PDB pdb(pdbFile,pdbList[i]);
        ProteinChain* pc = pdb.getFirstChain();
        string natSeq = pc->getSequence();
        vector<BackBoneSite*> bsList;
        proteinchain2BackboneSiteList(pc, bsList);

        DesignTemplate dt(bsList, &dp, s1Etable, s2Etable);
        dt.loadS1S2(s1Etable, s2Etable);
        dt.loadSingleResidueEnergy();
        dt.loadPairwiseEnergy();

        DesignMC mc(&dt);
        RotSequence unit(dt.resNum);
        mc.mcRun(&unit);
        string desSeq1 = unit.toAASequence();

        mc.mcRun(&unit);
        string desSeq2 = unit.toAASequence();

        mc.mcRun(&unit);
        string desSeq3 = unit.toAASequence();
        output << desSeq1 << endl;
        output << desSeq2 << endl;
        output << desSeq3 << endl;
        output.close();

    }


}



/*
 * Test40PDBAADesign.cpp
 *
 *  Created on: 2018Äê4ÔÂ10ÈÕ
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


    //double wtS2 = atof(args[1]);


    double totP1 = 0;
    double totP2 = 0;
    double s1TotP1 = 0;
    double s1TotP2 = 0;
    double s1TotRank = 0;
    int totRes = 0;


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

    double totRank = 0;

    for(int x=0;x<pdbList.size();x++){
        string pdbID = pdbList.at(x);
        string pdbFile = "/home/xiongpeng/sefTest/pdb/"+pdbID+".pdb";

        PDB pdb(pdbFile,pdbID);
        ProteinChain* pc = pdb.getFirstChain();

        string natSeq = pc->getSequence();

        vector<BackBoneSite*> bsList;
        proteinchain2BackboneSiteList(pc, bsList);


        AADesignTemplate adt2(bsList, &s1Etable, &s2Etable);

        /*
        AADesignTemplate adt2(bsList);
        string paFile = "/home/xiongpeng/abacus2.0/directPA/"+pdbID+".sa";
        string smFile = "/home/xiongpeng/abacus2.0/directSM/"+pdbID+".sm";


        adt2.loadS1S2(paFile,smFile);
*/


        string des = adt2.mcMinimize();
        SeqProfile prof1;
        adt2.getProfile(&prof1);

        vector<int> ranks = adt2.getNativeRankList();

        vector<int> s1Ranks = adt2.getNativeRankListS1();
        vector<double> s1PList = adt2.getNativePListS1();

        ResName rn;

        int len = natSeq.length();

        for(int i=0;i<len;i++){
            int aa = rn.sinToInt(natSeq[i]);
            double p1 = prof1.getP(i, aa);
            double p2 = 0;
            if(natSeq[i] == des[i])
                p2 = 1;
            totP1 += p1;
            totP2 += p2;
            s1TotP1 += s1PList[i];
            if(s1Ranks[i] == 1)
                s1TotP2 += 1.0;

            totRank += ranks[i];
            totRes++;
        }
    }

    double meanP1 = totP1/totRes;
    double meanP2 = totP2/totRes;

    double meanS1P1 = s1TotP1/totRes;
    double meanS1P2 = s1TotP2/totRes;
    double meanRank = totRank/totRes;

//    printf("wtLocal %4.2f wtLongRange: %4.2f meanRank: %5.3f\n", wtS2,0.8, meanRank);
    printf("s1 P1 P2: %7.5f %7.5f mean p1 p2: %7.5f %7.5f Rank: %6.4f\n",meanS1P1, meanS1P2, meanP1, meanP2, meanRank);

}

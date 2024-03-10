/*
 * TrainDesign.cpp
 *
 *  Created on: 2017Äê12ÔÂ27ÈÕ
 *      Author: notxp
 */


#include <iostream>
#include <vector>
#include <string>
#include <time.h>
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
#include "designseq/CmdArgs.h"

using namespace std;
using namespace NSPdesignseq;
using namespace NSPgeometry;

int main(int argc, char** args){
    /*
     * input is PDB structure
     */

    string pdbID = string(args[1]);
    string paraTag = string(args[2]);

    string pdbFile = "/export/home/s2982206/aspen/proteinDesign/trainingSet/pdb/"+pdbID+".pdb";
    string paraFile = "/export/home/s2982206/aspen/proteinDesign/trainingSet/para/" + paraTag;
    string output = "/export/home/s2982206/aspen/proteinDesign/trainingSet/designs/"+pdbID+"_"+paraTag+".pdb";

    DesignParameters dp(paraFile);

    PDB pdb(pdbFile, pdbID);


    vector<BackBoneSite*> bsList;
    pdb2BackboneSiteList(&pdb, bsList);

    S1EnergyTable s1Etable;
    S2MatrixFinder s2Etable;

    DesignTemplate* dt = new DesignTemplate(bsList, &dp, s1Etable, s2Etable,1.0);

    cout << "start load S1S2" << endl;

    dt->loadS1S2(s1Etable, s2Etable);
    cout << "start load single Residue Energy" << endl;
    dt->loadSingleResidueEnergy();
    cout << "start load pairwise energy" << endl;
    dt->loadPairwiseEnergy();


    cout << "start MC:" << endl;

    DesignMC mc(dt);
    RotSequence unit(dt->resNum);



    for(int i=1;i<=5;i++)
    {
        char x[10];
        char s[100];
        sprintf(x,"000%d",i);
        string xx(x);
        xx = xx.substr(xx.length()-3,3);
        string fileName = output;
        fileName = fileName.substr(0, output.length()-4) + "-" + xx + ".pdb";
        double score = mc.mcRun(&unit);
        mc.printPDB(&unit, fileName);
    }
    delete dt;

}


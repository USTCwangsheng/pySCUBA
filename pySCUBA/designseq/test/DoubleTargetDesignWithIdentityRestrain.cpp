/*
 * DoubleTargetDesignWithIdentityRestrain.cpp
 *
 *  Created on: 2018��3��15��
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

    CmdArgs cmd(argc, args);

    if(!cmd.specifiedOption("-in1") || !cmd.specifiedOption("-in2") || !cmd.specifiedOption("-out") || !cmd.specifiedOption("-idt") || !cmd.specifiedOption("-idtF")){
        cout << "Usage:\nDesignSeq -in1 $INPUTPDB1 -in2 $INPUTPDB2 -out $OUTPUTPDB -idt $IDENTITY -idtF $IDENTITYFACTOR" << endl;
        cout << "DesignSeq -in1 $INPUTPDB1 -in2 $INPUTPDB2 -out $OUTPUTPDB -n $SEQNUM -idt $IDENTITY -idtF $IDENTITYFACTOR" << endl;
        cout << "DesignSeq -in1 $INPUTPDB1 -in2 $INPUTPDB2 -out $OUTPUTPDB -resFile $RESFILE -idt $IDENTITY -idtF $IDENTITYFACTOR" << endl;
    }

    int n = 1;
    if(cmd.specifiedOption("-n")){
        n = atoi(cmd.getValue("-n").c_str());
    }

    double idt = atof(cmd.getValue("-idt").c_str());
    double idtF = atof(cmd.getValue("-idtF").c_str());

    cout << "start!" << endl;
    clock_t start = clock();
    string s1 = cmd.getValue("-in1");
    string s2 = cmd.getValue("-in2");
    string pdbID = "unk";
    string output = cmd.getValue("-out");

    if(output.length() < 4 || output.substr(output.length()-4, 4) != ".pdb")
        output = output + ".pdb";


    PDB pdb1(s1, pdbID);
    ProteinChain* pc1 = pdb1.getFirstChain();
    string natSeq1 = pc1->getSequence();

    vector<BackBoneSite*> bsList1;
    proteinchain2BackboneSiteList(pc1, bsList1);



    PDB pdb2(s2, pdbID);
    ProteinChain* pc2 = pdb2.getFirstChain();
    string natSeq2 = pc2->getSequence();

    vector<BackBoneSite*> bsList2;
    proteinchain2BackboneSiteList(pc2, bsList2);


    cout << "read model1 and model2" << endl;

    string path = NSPdataio::datapath();
    string refFile = path + "ref/hecRefDoubleTarget";
    DesignParameters dp;
    dp.readHECRef(refFile);


    S1EnergyTable s1Etable;
    S2MatrixFinder s2Etable;

    DesignTemplate* dt1 = new DesignTemplate(bsList1, &dp, s1Etable, s2Etable);
    DesignTemplate* dt2 = new DesignTemplate(bsList2, &dp, s1Etable, s2Etable);

    DoubleTargetDesign dmc(dt1, dt2, 1.0, 1.0);
    cout << "generate double target design template: " << endl;

    dmc.loadEnergyTable(s1Etable, s2Etable);

    cout << "energy table loaded" << endl;


    RotSequence unit1(dt1->resNum);
    RotSequence unit2(dt2->resNum);
    RotSequence unit3(dt1->resNum);
    RotSequence unit4(dt2->resNum);

    cout << "designing: " << endl;
    for(int i=1;i<=n;i++)
    {
        char x[10];
        sprintf(x,"000%d",i);
        string xx(x);
        xx = xx.substr(xx.length()-3,3);
        string fileNameA;
        string fileNameB;
        string fileNameC;
        string fileNameD;
        if(n == 1){
            fileNameA = output.substr(0, output.length()-4) + "-A" + ".pdb";
            fileNameB = output.substr(0, output.length()-4) + "-B" + ".pdb";
            fileNameC = output.substr(0, output.length()-4) + "-C" + ".pdb";
            fileNameD = output.substr(0, output.length()-4) + "-D" + ".pdb";
        }
        else {
            fileNameA = output.substr(0, output.length()-4) + "-A-" + xx + ".pdb";
            fileNameB = output.substr(0, output.length()-4) + "-B-" + xx + ".pdb";
            fileNameC = output.substr(0, output.length()-4) + "-C-" + xx + ".pdb";
            fileNameD = output.substr(0, output.length()-4) + "-D-" + xx + ".pdb";
        }
        dmc.threeTargetDesignMC(&unit1, &unit2, &unit3, &unit4, idt, idtF);
        dmc.printPDB(&unit1, fileNameA, dt1);
        dmc.printPDB(&unit2, fileNameB, dt2);
        dmc.printPDB(&unit3, fileNameC, dt1);
        dmc.printPDB(&unit4, fileNameD, dt2);
    }




    delete dt1;
    delete dt2;



    clock_t end = clock();
    cout << "time: " << (float)(end-start)/CLOCKS_PER_SEC << endl;

}

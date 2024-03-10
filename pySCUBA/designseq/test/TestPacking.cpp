/*
 * TestPacking.cpp
 *
 *  Created on: 2018Äê4ÔÂ24ÈÕ
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
#include "designseq/PackingTemplate.h"

using namespace std;
using namespace NSPdesignseq;
using namespace NSPgeometry;

int main(int argc, char** args){
    /*
     * input is PDB structure
     */
    CmdArgs cmd(argc, args);

    if(!cmd.specifiedOption("-in") || !cmd.specifiedOption("-out")){
        cout << "Usage:\nPacking -in $INPUTPDB -out $OUTPUTPDB" << endl;
        cout << "Packing -in $INPUTPDB -seq $SEQFILE -out $OUTPUTPDB" << endl;
        cout << "Packing -in $INPUTPDB -para $PARAFILE -out $OUTPUTPDB" << endl;
    }


    string s = cmd.getValue("-in");
    string pdbID = "unk";
    string output = cmd.getValue("-out");

//    string paraFile = "/home/xiongpeng/packing/para/" + cmd.getValue("-p");

    if(output.length() < 4 || output.substr(output.length()-4, 4) != ".pdb")
        output = output + ".pdb";

    PDB pdb(s, pdbID);
    vector<Residue*> resList = pdb.getResList();


    ResName rn;
    if(cmd.specifiedOption("-seq")){
        string seqFile = cmd.getValue("-seq");
        ifstream input;
        string seq = "";
            input.open(seqFile.c_str(),ios::in);
            if (! input.is_open())
            {
                cout << "fail to open file " << seqFile << endl;
                exit (1);
            }
            string ss;
            while(getline(input,ss))
            {
                if(ss[0]=='>') continue;
                seq = seq + ss;
            }

        if(seq.length() != resList.size()){
            cout << "sequence length: " << seq.length() << " residue number: " << resList.size() << " not match" << endl;
        }
        for(int i=0;i<seq.length();i++){
            resList[i]->triName = rn.sinToTri(seq[i]);
        }
    }

    DesignParameters* dp;
    if(cmd.specifiedOption("-para")) {
        string paraFile = cmd.getValue("-para");
        dp = new DesignParameters(paraFile);
    }
    else
        dp = new DesignParameters(2);

    PackingTemplate pt(&pdb,dp);
    //pt.printDetail();

    DesignMC mc(&pt);
    RotSequence unit(pt.resNum);


    mc.packing2(&unit);
    double ene = unit.totEnergy(&pt);
    cout << "ene: " << ene << endl;
    mc.printPDB(&unit,output);

    delete dp;
}


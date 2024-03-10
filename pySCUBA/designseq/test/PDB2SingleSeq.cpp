/*
 * PDB2SingleSeq.cpp
 *
 *  Created on: Mar 21, 2019
 *      Author: hbfrank
 */

#include <iostream>
#include <vector>
#include <string>
#include "stdio.h"
#include "designseq/ProteinRep.h"


using namespace std;
using namespace NSPdesignseq;

int main(int argc, char** args){
    if (argc < 2) {
        cerr << "[Error] Missing PDB file!" << endl;
    }

    string pdbFile = string(args[1]);
    PDB* pdb = new PDB(pdbFile, "UNK");
    vector<ProteinChain*> pcList = pdb->getChains();
    for(int i=0;i<pcList.size();i++) {
        ProteinChain* pc = pcList[i];
        char c = pc->getChainID();
        string seq = pc->getSequence();
        cout << seq;
    }
    cout << endl;

    delete pdb;
}




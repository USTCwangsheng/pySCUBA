/*
 * PDB2Seqs.cpp
 *
 *  Created on: Oct 3, 2018
 *      Author: s2982206
 */

#include <iostream>
#include <vector>
#include <string>
#include "stdio.h"
#include "designseq/ProteinRep.h"


using namespace std;
using namespace NSPdesignseq;


int main(int argc, char** args){



    string pdbFile = string(args[1]);
    PDB* pdb = new PDB(pdbFile, "UNK");
    if(argc == 2) {
        vector<ProteinChain*> pcList = pdb->getChains();
        for(int i=0;i<pcList.size();i++) {
            ProteinChain* pc = pcList[i];
            char c = pc->getChainID();
            string seq = pc->getSequence();
            cout << c << " " << seq << endl;
        }
    }
    if(argc == 3) {
        char chainID = args[2][0];
        ProteinChain* pc = pdb->getChain(chainID);
        if(pc == NULL) {
            cout << "invalid chainID: " << args[2] << endl;
            exit(0);
        }
        string seq = pc->getSequence();
        cout << seq << endl;
    }
    delete pdb;
}




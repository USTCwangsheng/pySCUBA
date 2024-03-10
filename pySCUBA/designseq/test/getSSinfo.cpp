/*
 * getSSinfo.cpp
 *
 *  Created on: Apr 6, 2019
 *      Author: xuyang
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
#include "designseq/ScoringTemplate.h"

using namespace NSPdesignseq;
using namespace NSPgeometry;

int main(int argc, char**argv) {
    if(argc==1) {
        std::cout <<"Usage:" <<std::endl;
        std::cout <<"\t${pdbfile} ${outfile}" <<std::endl;
        return 0;
    }

    std::string pdbfile(argv[1]);
    std::string outfile(argv[2]);

    string pdbID = "unk";
    PDB pdb(pdbfile, pdbID);

    vector<ProteinChain*> pcs=pdb.getChains();
    vector<string> ssList;

    for(auto &pc:pcs) {
        StructureInfo* strInfo = new StructureInfo(pc);
        SasaPSD sasaPSD;
        strInfo->updateTorsion();
        strInfo->updateSecondaryStructure();
        string ssList1;
        int resNum = pc->getChainLength();
        for(int i=0;i<resNum;i++){
            char ss = strInfo->getSS(i);
            ssList1.push_back(ss);
        }
        ssList.push_back(ssList1);
        delete strInfo;
    }
    ofstream ofs(outfile);
    for(auto &s:ssList) ofs <<s <<std::endl;
    ofs.close();
}


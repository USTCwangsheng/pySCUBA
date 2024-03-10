/*
 * getBBHBinfo.cpp
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
    std::string pdbID = "unk";
    PDB pdb(pdbfile, pdbID);
    vector<ProteinChain*> pcs=pdb.getChains();

    std::ofstream ofs(argv[2]);
    for(int i=0;i<pcs.size();i++) {
        StructureInfo* strInfo = new StructureInfo(pcs[i]);
        strInfo->updateBBHBonds();
        std::vector<BackboneHBond*> bhbs = strInfo->gethblist();
        for(auto &h:bhbs) {
            ofs <<i <<std::endl;
            ofs <<h->donorseq <<std::endl;
            ofs <<i <<std::endl;
            ofs <<h->accepterseq <<std::endl;
        }
        delete strInfo;
    }
    ofs.close();
}


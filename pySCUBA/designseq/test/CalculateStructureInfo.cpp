/*
 * CalculateStructureInfo.cpp
 *
 *  Created on: 2018Äê3ÔÂ19ÈÕ
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

using namespace std;
using namespace NSPdesignseq;
using namespace NSPgeometry;

int main(int argc, char** args){

    /*
     * input is PDB structure
     */
    if(argc != 2)
    {
        cout << "structInfo $PDBFile" << endl;
    }
    string s(args[1]);
    string pdbID = "unk";
    PDB pdb(s, pdbID);
    ProteinChain* pc = pdb.getFirstChain();

    /*
     * calculate backbone torsion angle, secondary structure, and SASA
     */
    StructureInfo strInfo(&pdb);

    SasaPSD sasaPSD;
    strInfo.updateTorsion();
    strInfo.updateSecondaryStructure();
    strInfo.updateSAI(&sasaPSD);

    int resNum = strInfo.getResNum();

    cout << "C RES SEQ SS   PHI     PSI     OMG    SAI" << endl;
    for(int i=0;i<resNum;i++){
        Residue* res = strInfo.getResidue(i);
        string resID = res->resID;
        char chainID = res->chainID;
        char ss = strInfo.getSS(i);
        float phi = strInfo.getPhi(i);
        float psi = strInfo.getPsi(i);
        float omg = strInfo.getOmg(i);
        float sai = strInfo.getSai(i);
        string triName = res->triName;
        printf("%c %3s %s  %c %7.2f %7.2f %7.2f %5.2f\n", chainID, resID.c_str(), triName.c_str(), ss, phi, psi, omg, sai);
    }
}


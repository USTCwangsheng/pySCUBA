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
#include "designseq/PackingTemplate.h"

using namespace std;
using namespace NSPdesignseq;
using namespace NSPgeometry;

int main(int argc, char** args){

    /*
     * input is PDB structure
     */
    cout << "start" << endl;

    if(argc != 2)
    {
        cout << "testScoring $PDBFile" << endl;
    }
    string s(args[1]);
    string pdbID = "unk";
    PDB pdb(s, pdbID);
    ProteinChain* pc = pdb.getFirstChain();

    /*
     * calculate backbone torsion angle, secondary structure, and SASA
     */
    StructureInfo* strInfo = new StructureInfo(pc);

    SasaPSD sasaPSD;
    strInfo->updateTorsion();
    strInfo->updateSecondaryStructure();
    strInfo->updateSAI(&sasaPSD);


    vector<float> saiList;
    vector<float> s1List;
    vector<float> s2List;
    vector<float> packList;
    vector<char> ssList;
    vector<float> rotEList;

    double totS1 = 0;
    /*
     * print local structure information
     */
//    cout << "local structure information:" << endl;
    int resNum = pc->getChainLength();
    for(int i=0;i<resNum;i++){
        Residue* res = pc->getResList().at(i);
        string resID = res->resID;
        char ss = strInfo->getSS(i);
        float phi = strInfo->getPhi(i);
        float psi = strInfo->getPsi(i);
        float omg = strInfo->getOmg(i);
        float sai = strInfo->getSai(i);
        string triName = res->triName;
//        printf("%3s %s %c %7.2f %7.2f %7.2f %5.2f\n", resID.c_str(), triName.c_str(), ss, phi, psi, omg, sai);
        ssList.push_back(ss);
        saiList.push_back(sai);
        s1List.push_back(0);
        s2List.push_back(0);
        packList.push_back(0);
        rotEList.push_back(0);
    }

    /*
     * convert ProteinChain to BackboneSite list
     * in this step, the torsion angles, secondary structure, and SASA will be calculated
     */

    vector<BackBoneSite*> bsList;
    proteinchain2BackboneSiteList(pc, bsList);

    /*
     * calculate S1 score
     */

    S1EnergyTable s1Etable;
    ResName rn;
    cout << "S1 score: " << endl;

    for(int i=0;i<resNum;i++){
        Residue* res = pc->getResList().at(i);
        string resID = res->resID;
        string triName = res->triName;

        AAProbabilityArray pa;
        s1Etable.getS1(*bsList.at(i), &pa);
        int intName = rn.triToInt(triName);
        float s1 = pa.getScore(intName);
        totS1 += s1;
        s1List[i] = s1;
    }

    //    cout << "S2 score: " << endl;


        double totS2 = 0;

        /*
         * update residue pairs
         */
        S2MatrixFinder s2Etable;
        for(int i=0;i<resNum;i++){
            BackBoneSite* bsA = (bsList.at(i));
            string triNameA = pc->getResList().at(i)->triName;
            int intNameA = rn.triToInt(triNameA);
            string resIDA = pc->getResList().at(i)->resID;

            for(int j=i+1;j<resNum;j++){
                BackBoneSite* bsB = (bsList.at(j));
                string triNameB = pc->getResList().at(j)->triName;
                int intNameB = rn.triToInt(triNameB);
                string resIDB = pc->getResList().at(j)->resID;
                float cbDistance = bsA->cbcrd().distance(bsB->cbcrd());
                /*
                 * residue pair is defined by CB atom distance smaller than 8.0
                 */
                if(cbDistance > 8.0) continue;
                BackboneSitesPair pair((bsList.at(i)), (bsList.at(j)));
                AAScoreMatrix sm;
                s2Etable.getSM(&pair, &sm);
                int sep = j-i;
                float s2 = sm.getValue(intNameA, intNameB);


                if(sep < 5 || (ssList[i] == 'E' && ssList[j] == 'E'))
                    s2 = s2* 0.4;
                else
                    s2 = s2 * 0.7;
                totS2 += s2;

                s2List[i] += 0.5*s2;
                s2List[j] += 0.5*s2;
            }
        }

        /*
         * calculate atomic score
         */

        double totVdwEnergy = 0;

        AtomicEnergyCalcular ec;
        DesignParameters dp;
        AtomLib atLib;
        vector<Rotamer*> rotList;
        for(int i=0;i<resNum;i++){
            Residue* res = pc->getResList().at(i);
            rotList.push_back(res->natRotamer(&atLib));
        }

        for(int i=0;i<resNum;i++){
            string resIDA = pc->getResList().at(i)->resID;
            BackBoneSite* bsA = (bsList.at(i));
            Rotamer* rotA = rotList.at(i);
            string triNameA = rotA->triName;
            if(rotA->triName == "XXX") continue;

            Conformer* confA = new Conformer(rotA, bsA, &atLib);
            double saiA = strInfo->getSai(i);

            for(int j=i+2;j<resNum;j++){
                string resIDB = pc->getResList().at(j)->resID;

                double saiB = strInfo->getSai(j);
                float vdwWT = dp.getVdwWeight(0.5*saiA + 0.5*saiB);
                BackBoneSite* bsB = (bsList.at(j));
                Rotamer* rotB = rotList.at(j);
                string triNameB = rotB->triName;
                if(triNameB == "XXX") continue;
                Conformer* confB = new Conformer(rotB, bsB, &atLib);

                int seqSep = j-i;
                if(seqSep > 5) seqSep = 5;
                float e = pairEnergyPack(confA, confB, seqSep, vdwWT, &ec);

                packList[i] += 0.5*e;
                packList[j] += 0.5*e;

                totVdwEnergy += e;
            }
        }

        /*
         * calculate rotamer energy
         */

        double totRotEnergy = 0;

        //RotamerLib rotLibA("A1");
        //RotamerLib rotLibB("B1");


        //PhipsiLib ppLib;

        /*
        for(int i=0;i<resNum;i++){
            double phi = bsList[i]->phi();
            double psi = bsList[i]->psi();
            Phipsi pp(phi,psi);

            int index = ppLib.phipsiToIndex(&pp);
            double e = 0;
            if(pp.regionAB() == 'A')
                e = rotLibA.getNearestRotamerEnergy(rotList[i], index);
            else
                e = rotLibB.getNearestRotamerEnergy(rotList[i], index);

            totRotEnergy += e;

            rotEList[i] = e;
        }

        cout << "resID AA ss  SAI   S1    S2    VDW    EROT" << endl;
        */
}

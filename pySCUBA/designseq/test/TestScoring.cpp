/*
 * TestScoring.cpp
 *
 *  Created on: 2017��12��25��
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
#include "designseq/ScoringTemplate.h"


using namespace std;
using namespace NSPdesignseq;
using namespace NSPgeometry;

int main(int argc, char** args){

    /*
     * input is PDB structure
     */

    if(argc == 1)
    {
        std::cout <<"Usage:" <<std::endl;
        cout << "\t${PDBFile} es1=${s1_outfile}(optional) es2=${s2_outfile}(optional) "
                << "pck=${packing_outfile}(optional) rot=${rotamer_outfile}(optional) ref=${ref_outfile}(optional)" << endl;
        return 1;
    }
    std::map<std::string,std::string> fs;
    for(int i=2;i<argc;i++) {
        std::string f(args[i]);
        fs.insert({f.substr(0,3),f.substr(4)});
    }
    string s(args[1]);
    string pdbID = "unk";
    PDB pdb(s, pdbID);

    //string parafile(args[2]);
    //DesignParameters dp(parafile);
    DesignParameters dp;

    S1EnergyTable s1Etable;
    S2MatrixFinder s2Etable;
    ScoringTemplate* st = new ScoringTemplate(&pdb, &dp, s1Etable, s2Etable);


    /*
     * calculate backbone torsion angle, secondary structure, and SASA
     */
    StructureInfo* strInfo = new StructureInfo(&pdb);


    SasaPSD sasaPSD;
    strInfo->updateTorsion();
    strInfo->updateSecondaryStructure();
    strInfo->updateSAI(&sasaPSD);


    vector<float> saiList;
    vector<float> s1List;
    vector<float> s2List;
    std::map<std::vector<int>,float> s2s;
    vector<float> packList;
    std::map<std::vector<int>,float> pks;
    vector<float> refList;
    vector<char> ssList;
    vector<float> rotEList;



    /*
     * print local structure information
     */
//    cout << "local structure information:" << endl;
    int resNum = strInfo->getResNum();
    for(int i=0;i<resNum;i++){
        Residue* res = strInfo->getResidue(i);
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
        refList.push_back(0);
    }



    /*
     * calculate S1 score
     */

    ResName rn;

    float totS1 = 0;
    float totRef = 0;
    float totRot = 0;

    for(int i=0;i<resNum;i++){
        Residue* res = strInfo->getResidue(i);
        string resID = res->resID;
        string triName = res->triName;
        float s1 = st->s1EAList[i]->getEnergy(0);
        float ref = st->refEAList[i]->getEnergy(0);
        float rot = st->rotEAList[i]->getEnergy(0);
        totS1 += s1;
        totRef += ref;
        totRot += rot;
        s1List[i] = s1;
        refList[i] = ref;
        rotEList[i] = rot;
    }


    double totS2 = 0;
    double totVdwEnergy = 0;

    for(int i=0;i<st->s2EMList.size();i++){
        BackboneSitesPair* pair = st->bsPairs[i];
        EnergyMatrix* em = st->s2EMList[i];
        int posA = pair->siteA->resseq;
        int posB = pair->siteB->resseq;
        float s2 = em->getEnergy(0,0);
        totS2 += s2;
        s2List[posA] += 0.5*s2;
        s2List[posB] += 0.5*s2;
        s2s.insert({{posA,posB},s2});
    }

    for(int i=0;i<st->packEMList.size();i++){
        BackboneSitesPair* pair = st->bsPairs[i];
        EnergyMatrix* em = st->packEMList[i];
        int posA = pair->siteA->resseq;
        int posB = pair->siteB->resseq;
        float e = em->getEnergy(0,0);
        packList[posA] += 0.5*e;
        packList[posB] += 0.5*e;
        totVdwEnergy += e;
        pks.insert({{posA,posB},e});
    }


    cout << "resID AA ss  SAI   S1    S2    VDW    EROT    REF" << endl;

    for(int i=0;i<resNum;i++){
        Residue* res = strInfo->getResidue(i);
        string resID = res->resID;
        string type = res->triName;
        char ss = ssList[i];
        float sai = saiList[i];
        float s1 = s1List[i];
        float s2 = s2List[i];
        float vdw = packList[i];
        float rot = rotEList[i];
        float ref = refList[i];

        printf("%-4s %3s %c %4.2f %5.2f %6.2f %6.2f %6.2f %6.2f\n",resID.c_str(),type.c_str(),ss,sai,s1,s2,vdw,rot,ref);
    }

    cout << endl;
    double totEnergy = totS1 + totS2 + totVdwEnergy + totRef + totRot;
    printf("TOTAL S1: %6.2f S2: %6.2f PACKING: %6.2f ROT: %6.2f REF: %6.2f TOT: %6.2f\n", totS1, totS2, totVdwEnergy, totRot, totRef, totEnergy);

    for(auto &f:fs) {
        if(f.first=="es1") {
            std::ofstream ofs(f.second);
            for(auto &e:s1List) ofs <<e <<std::endl;
            ofs.close();
        } else if(f.first=="es2") {
            std::ofstream ofs(f.second);
            for(auto &p:s2s) {
                ofs <<p.first[0] <<std::endl;
                ofs <<p.first[1] <<std::endl;
                ofs <<p.second <<std::endl;
            }
            ofs.close();
        } else if(f.first=="pck") {
            std::ofstream ofs(f.second);
            for(auto &p:pks) {
                ofs <<p.first[0] <<std::endl;
                ofs <<p.first[1] <<std::endl;
                ofs <<p.second <<std::endl;
            }
            ofs.close();
        } else if(f.first=="rot") {
            std::ofstream ofs(f.second);
            for(auto &e:rotEList) ofs <<e <<std::endl;
            ofs.close();
        } else if(f.first=="ref") {
            std::ofstream ofs(f.second);
            for(auto &e:refList) ofs <<e <<std::endl;
            ofs.close();
        }
    }
}



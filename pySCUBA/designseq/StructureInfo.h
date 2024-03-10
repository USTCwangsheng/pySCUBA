/*
 * StructureInfo.h
 *
 *  Created on: 2017年11月1日
 *      Author: notxp
 */

#ifndef DESIGNSEQ_STRUCTUREINFO_H_
#define DESIGNSEQ_STRUCTUREINFO_H_

#include <vector>
#include <string>
#include <map>
#include "designseq/ProteinRep.h"
#include "designseq/Angles.h"
#include "designseq/SasaPSD.h"

namespace NSPdesignseq {

using namespace NSPproteinrep;

class BackboneHBond{
public:
    char donerChainID;
    char acceptorChainID;
    string donerResID;
    string acceptorResID;
    int seqSeparation;
    XYZ hDoner;
    XYZ hydrogen;
    XYZ hAcceptor;
    float distance;
    float angle;

    BackboneHBond(XYZ& N, XYZ& H, XYZ& O, Residue* resA, Residue* resB);
    string toString();


    /*
     * added by xuyang, 2019.4.6
     */
    int donorseq;
    int accepterseq;



    virtual ~BackboneHBond(){}
};


class StructureInfo {
private:
    int resNum;
    vector<Residue*> resList;
    vector<XYZ*> NList;
    vector<XYZ*> CAList;
    vector<XYZ*> CList;
    char* ssSeq;
    float* phiList;
    float* psiList;
    float* omgList;
    map<string,int> chainIDResIDToSeqID;
    vector<BackboneHBond*> hbList; //need delete
    vector<float> saiList;
    vector<bool> backboneHBonded;



public:
    StructureInfo(PDB* protein);
    StructureInfo(ProteinChain* chain);
    StructureInfo(vector<Residue*>& residueList);

    void updateTorsion();
    void updateBBHBonds();
    void updateSecondaryStructure();
    void updateSAI(SasaPSD* rsp);

    string getSecondaryStructureSeq();
    Phipsi getPhipsi(Residue* res);
    float getSai(Residue* res);

    int getResNum() {return this->resNum;}
    Residue* getResidue(int seqID){
        if(seqID < 0 || seqID >= resNum){
            cout << "out of index: " << seqID << endl;
            exit(1);
        }
        return resList.at(seqID);
    }

    string getResID(int seqID) {
        if(seqID < 0 || seqID >= resNum){
            cout << "out of index: " << seqID << endl;
            exit(1);
        }
        return resList.at(seqID)->resID;
    }


    float getPhi(int seqID) {
        if(seqID < 0 || seqID >= resNum){
            cout << "out of index: " << seqID << endl;
            exit(1);
        }
        return phiList[seqID];
    }
    float getPsi(int seqID) {
        if(seqID < 0 || seqID >= resNum){
            cout << "out of index: " << seqID << endl;
            exit(1);
        }
        return psiList[seqID];
    }
    float getOmg(int seqID) {
        if(seqID < 0 || seqID >= resNum){
            cout << "out of index: " << seqID << endl;
            exit(1);
        }
        return omgList[seqID];
    }
    char getSS(int seqID) {
        if(seqID < 0 || seqID >= resNum){
            cout << "out of index: " << seqID << endl;
            exit(1);
        }
        char c = ssSeq[seqID];
        if(c == 'H' || c == 'G')
            return 'H';
        else if(c == 'C')
            return 'C';
        else
            return 'E';
    }

    float getSai(int seqID){
        if(seqID < 0 || seqID >= resNum){
            cout << "out of index: " << seqID << endl;
            exit(1);
        }
        return saiList.at(seqID);
    }

    int getFlexibleIndex(int seqID){
        if(this->backboneHBonded.size() != resNum){
            cout << "backbone hbond information insufficient" << endl;
            abort();
        }

        if(seqID < 0 || seqID >= resNum){
            cout << "out of index: " << seqID << endl;
            exit(1);
        }
        int fi = 0;
        if(seqID-2 >=0 && this->backboneHBonded[seqID-2])
            fi += 1;
        if(seqID-1 >=0 && this->backboneHBonded[seqID-1])
            fi += 2;
        if(this->backboneHBonded[seqID])
            fi += 4;
        if(seqID+1 < resNum && this->backboneHBonded[seqID+1])
            fi += 2;
        if(seqID+2 < resNum && this->backboneHBonded[seqID+2])
            fi += 1;
        return fi;
    }

    vector<float> getSaiList(){
        return saiList;
    }


    /*
     * added by xuyang, 2019.4.6
     */
    vector<BackboneHBond*> gethblist() {return hbList;}

    virtual ~StructureInfo();
};

void backboneSiteListUpdateSasa(vector<BackBoneSite*>& bsList);

void proteinchain2BackboneSiteList(ProteinChain* pc, vector<BackBoneSite*>& bsList);

void pdb2BackboneSiteList(PDB* pdb, vector<BackBoneSite*>& bsList);

void resList2BackboneSiteList(vector<Residue*>& resList, vector<NSPproteinrep::BackBoneSite>& bsList);




/*
 * added by xuyang, 2019.3.31
 */

void pdb2BackboneSiteList(PDB* pdb, vector<BackBoneSite*>& bsList, std::string &ssseq);


} /* namespace NSPdesignseq */

#endif /* DESIGNSEQ_STRUCTUREINFO_H_ */

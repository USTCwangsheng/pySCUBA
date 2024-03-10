/*
 * SingleSitePackingTemplate.cpp
 *
 *  Created on: Sep 21, 2018
 *      Author: s2982206
 */

#include "designseq/SingleSitePackingTemplate.h"

namespace NSPdesignseq {

SingleSitePackingTemplate::SingleSitePackingTemplate(ProteinChain* pc, DesignParameters* dp, char mutChainID, string mutResID, string mutType) {
    // TODO Auto-generated constructor stub

    vector<Residue*> tmpResList = pc->getResList();
    for(unsigned int i=0;i<tmpResList.size();i++){
        Residue* res = tmpResList.at(i);
        if(res->triName == "MSE") {
            res->triName = "MET";
            Atom* a = res->getAtom("SE");
            if(a != NULL)
                a->name = "SD";
        }

        if(res->hasThreeCoreAtoms() && res->hasAtom("O")){
            this->resList.push_back(res);
            Rotamer* rot = NULL;
            if(res->sidechainComplete(&atLib))
                rot = res->natRotamer(&atLib);
            this->natRotList.push_back(rot);
        }
    }


    proteinchain2BackboneSiteList(pc, this->bsList);


    this->resNum = bsList.size();

    this->aec = AtomicEnergyCalcular(dp);
    char s[20];
    char ch = 'A';


    for(int i=0;i<resNum;i++){
        BackBoneSite* bs = bsList.at(i);
        string resID = bs->strresid;
        ch = bs->chainid;
        sprintf(s,"%c%s",ch,resID.c_str());
        string key = string(s);
        chainIDResIDToIndex[key] = i;
        RotamerGroup* gp = new RotamerGroup();
        rotGroups.push_back(gp);
    }

    sprintf(s,"%c%s",mutChainID,mutResID.c_str());
    string mutkey = string(s);

    if(chainIDResIDToIndex.count(mutkey) == 0) {
        cerr << "invalid residue: " << mutkey << endl;
    }

    this->mutPosition = chainIDResIDToIndex[mutkey];
    this->mutType = mutType;
    this->dp = dp;
    this->rotLibA = new RotamerLib("A2");
    this->rotLibB = new RotamerLib("B2");


    nativeAARotamerGroup();
    updateResPairs();
    nonContactRotamer();
    loadSingleResidueEnergy();
    loadPairwiseEnergy();
}

void SingleSitePackingTemplate::nativeAARotamerGroup() {
    PhipsiLib ppLib;
    for(int pos=0;pos<this->resNum;pos++){
        float phi = bsList.at(pos)->phi();
        float psi = bsList.at(pos)->psi();

        Phipsi pp(phi, psi);
        int ppType = ppLib.phipsiToIndex(&pp);
        RotamerLib* rotLib;
        if(pp.regionAB() == 'A')
            rotLib = rotLibA;
        else
            rotLib = rotLibB;
        RotamerGroup* gp = rotGroups[pos];
        string triName = bsList[pos]->resname;
        if(pos == mutPosition) {
            triName = mutType;
        }
        if(triName == "MSE")
            triName = "MET";
        RotamerGroup* gp0 = rotLib->getAAGroup(triName);
        for(Rotamer* rot : gp0->rotList) {
            gp->addRotamer(rot);
        }
    }
}

void SingleSitePackingTemplate::nonContactRotamer(){
    vector<int>* mutNeighbors = neighborResidues[mutPosition];
    set<int> neighbors;
    for(int i=0;i<mutNeighbors->size();i++){
        neighbors.insert(mutNeighbors->at(i));
    }

    for(int pos=0;pos<this->resNum;pos++) {
        if(pos == mutPosition) continue;
        if(neighbors.find(pos) != neighbors.end()) continue;
        if(this->natRotList[pos] == NULL) continue;
        this->fixedPos.insert(pos);
        RotamerGroup* gp = rotGroups[pos];
        gp->clear();
        gp->addRotamer(this->natRotList[pos]);
    }

}

void SingleSitePackingTemplate::updateResPairs(){
    /*
     * update mutPosition involved bsPairs
     */
    this->bsPairs.clear();
    for(int i=0;i<this->resNum;i++){
        BackBoneSite* bsA = bsList[i];
        XYZ cbA = bsA->cbcrd();
        for(int j=i+1;j<this->resNum;j++){
            BackBoneSite* bsB = bsList[j];
            XYZ cbB = bsB->cbcrd();
            double d = cbA.distance(cbB);
            if(d < 10.5)
                bsPairs.push_back(new BackboneSitesPair(bsA, bsB));
        }
    }

    /*
     * update involvedPairs and neighbor residues
     */

    for(int i=0;i<involvedPairs.size();i++){
        delete involvedPairs.at(i);
    }

    for(int i=0;i<neighborResidues.size();i++){
        delete neighborResidues.at(i);
    }

    this->involvedPairs.clear();
    this->neighborResidues.clear();


    for(int i=0;i<resNum;i++){
        this->involvedPairs.push_back(new vector<int>());
        this->neighborResidues.push_back(new vector<int>());
    }

    for(int i=0;i<bsPairs.size();i++){
        BackboneSitesPair* pair = bsPairs[i];
        int idA = pair->siteA->resseq;
        int idB = pair->siteB->resseq;
        involvedPairs.at(idA)->push_back(i);
        involvedPairs.at(idB)->push_back(i);
        neighborResidues.at(idA)->push_back(idB);
        neighborResidues.at(idB)->push_back(idA);
    }
}

bool SingleSitePackingTemplate::contact(int posA, int posB, float cutoff){
    BackBoneSite* bsA = bsList[posA];
    BackBoneSite* bsB = bsList[posB];
    LocalFrame csA = getBackboneSiteLocalFrame(*bsA);
    LocalFrame csB = getBackboneSiteLocalFrame(*bsB);
    RotamerGroup* groupA = this->rotGroups[posA];
    RotamerGroup* groupB = this->rotGroups[posB];

    vector<XYZ> bbAList;
    vector<XYZ> bbBList;
    bbAList.push_back(bsA->ncrd());
    bbAList.push_back(bsA->cacrd());
    bbAList.push_back(bsA->ccrd());
    bbAList.push_back(bsA->ocrd());

    bbBList.push_back(bsB->ncrd());
    bbBList.push_back(bsB->cacrd());
    bbBList.push_back(bsB->ccrd());
    bbBList.push_back(bsB->ocrd());

    vector<vector<XYZ>> scAListList;
    vector<vector<XYZ>> scBListList;

    for(int i=0;i<groupA->rotNum;i++){
        Rotamer* rotA = groupA->rotList[i];
        scAListList.push_back(rotA->coordList);
    }

    for(int i=0;i<groupB->rotNum;i++){
        Rotamer* rotB = groupB->rotList[i];
        scBListList.push_back(rotB->coordList);
    }

    for(int i=0;i<groupA->rotNum;i++){
        for(XYZ a : scAListList[i]){
            for(XYZ b : bbBList){
                if(a.distance(b) < cutoff)
                    return true;
            }
        }
    }

    for(int i=0;i<groupB->rotNum;i++){
        for(XYZ a : scBListList[i]){
            for(XYZ b : bbAList){
                if(a.distance(b) < cutoff)
                    return true;
            }
        }
    }

    for(int i=0;i<groupA->rotNum;i++){
        for(int j=0;j<groupB->rotNum;j++){
            for(XYZ a : scAListList[i]){
                for(XYZ b : scBListList[j]){
                    if(a.distance(b) < cutoff)
                        return true;
                }
            }
        }
    }

    return false;
}

void SingleSitePackingTemplate::loadSingleResidueEnergy() {
    PhipsiLib ppLib;
    for(int i=0;i<resNum;i++) {
        RotamerGroup* group = rotGroups.at(i);
        int rotNum = group->rotNum;
        float phi = bsList.at(i)->phi();
        float psi = bsList.at(i)->psi();
        Phipsi pp(phi,psi);
        int ppIndex = ppLib.phipsiToIndex(&pp);
        char ppType = pp.regionAB();
        EnergyArray* ea = new EnergyArray(rotNum);
        if(rotNum == 1) {
            ea->setEnergy(0,0.0);
            this->eaList.push_back(ea);
            continue;
        }

        for(int j=0;j<rotNum;j++){
            Rotamer* rot = group->rotList.at(j);
            int aaType = rn.triToInt(rot->triName);
            float eRot;
            if(ppType == 'A')
                eRot = rotLibA->getRotamerEnergy(rot->rotName, ppIndex);
            else
                eRot = rotLibB->getRotamerEnergy(rot->rotName, ppIndex);

            ea->setEnergy(j, eRot);
        }
        this->eaList.push_back(ea);
    }
}

void SingleSitePackingTemplate::loadPairwiseEnergy() {
    vector<ConformerGroup*> cgList;
    for(int i=0;i<resNum;i++){
        cgList.push_back(new ConformerGroup(rotGroups.at(i), bsList.at(i), &atLib));
    }
    for(int i=0;i<this->bsPairs.size();i++){

        BackboneSitesPair* pair = this->bsPairs[i];
        int idA = pair->siteA->resseq;
        int idB = pair->siteB->resseq;
        int seqSep = abs(idA-idB);
        if(seqSep > 5) seqSep = 5;
        if(pair->siteA->chainid != pair->siteB->chainid) seqSep = 5;

        ConformerGroup* groupA = cgList[idA];
        ConformerGroup* groupB = cgList[idB];

        EnergyMatrix* em = new EnergyMatrix(groupA->rotNum, groupB->rotNum, idA, idB);

        if(this->fixedPos.find(idA) != this->fixedPos.end() && this->fixedPos.find(idB) != this->fixedPos.end()) {
            em->setEnergy(0, 0, 0.0);
            this->emList.push_back(em);
            continue;
        }

        for(int j=0;j<groupA->rotNum;j++){
            Conformer* confA = groupA->confList[j];
            int aaTypeA = rn.triToInt(confA->triName);
            for(int k=0;k<groupB->rotNum;k++){
                Conformer* confB = groupB->confList[k];
                int aaTypeB = rn.triToInt(confB->triName);
                float vdw = pairEnergyPack(confA, confB, seqSep, 1.0, &aec);
                em->setEnergy(j,k, vdw);
            }
        }
        this->emList.push_back(em);
    }

    for(int i=0;i<cgList.size();i++){
        delete cgList[i];
    }
}

void SingleSitePackingTemplate::printChoice() {
    for(int i=0;i<resNum;i++) {
        RotamerGroup* gp = this->rotGroups[i];
        int choiceNum = gp->rotNum;
        BackBoneSite* bs = bsList.at(i);
    //    cout << "pos: " << i << " choiceNum: " << choiceNum <<" " << bs->resseq<< endl;
    }

    int eaSize = this->eaList.size();
    int emSize = this->emList.size();
    int bsSize = this->bsPairs.size();
    int invSize = this->involvedPairs.size();
    int gpSize = this->rotGroups.size();

    printf("ea: %d em: %d bs: %d inv: %d gp: %d\n",eaSize, emSize, bsSize, invSize, gpSize);

}

SingleSitePackingTemplate::~SingleSitePackingTemplate() {

    delete rotLibA;
    delete rotLibB;


    for(int i=0;i<eaList.size();i++){
        delete eaList.at(i);
    }

    for(int i=0;i<rotGroups.size();i++){
        delete rotGroups.at(i);
    }

    for(int i=0;i<natRotList.size();i++) {
        delete natRotList.at(i);
    }

    for(int i=0;i<emList.size();i++){
        delete emList.at(i);
    }

    for(int i=0;i<involvedPairs.size();i++){
        delete involvedPairs.at(i);
    }

    for(int i=0;i<neighborResidues.size();i++){
        delete neighborResidues.at(i);
    }

    // TODO Auto-generated destructor stub
}

}

/*
 * ScoringTemplate.cpp
 *
 *  Created on: Dec 20, 2018
 *      Author: s2982206
 */

#include "designseq/ScoringTemplate.h"

namespace NSPdesignseq {

ScoringTemplate::ScoringTemplate(PDB* targetPDB, DesignParameters* dp, S1EnergyTable& s1Etable, S2MatrixFinder& s2Etable) {
    StructureInfo si(targetPDB);
    for(int i=0;i<si.getResNum();i++){
        this->resList.push_back(si.getResidue(i));
    }

    vector<BackBoneSite*> bsList;
    pdb2BackboneSiteList(targetPDB, bsList);
    backboneSiteListUpdateSasa(bsList);
    this->s1ET = &s1Etable;
    this->resNum = bsList.size();
    this->aec = new AtomicEnergyCalcular(dp);
    char s[20];
    char chainID = 'A';

    for(int i=0;i<bsList.size();i++){
        BackBoneSite* bs = bsList.at(i);
        this->bsList.push_back(bs);
        string resID = bs->strresid;
        chainID = bs->chainid;
        sprintf(s,"%c%s",chainID,resID.c_str());
        string key = string(s);
        chainIDResIDToIndex[key] = i;
        RotamerGroup* gp = new RotamerGroup();
        rotGroups.push_back(gp);
    }
    this->dp = dp;

    this->rotLibA = new RotamerLib("A1");
    this->rotLibB = new RotamerLib("B1");

    setRotamerGroups();
    updateResPairs();
    loadS1S2(s1Etable, s2Etable);
    loadSingleResidueEnergy();
    loadPairwiseEnergy();
}

void ScoringTemplate::loadSingleResidueEnergy() {
    PhipsiLib ppLib;
    for(int i=0;i<resNum;i++){
        RotamerGroup* group = rotGroups.at(i);
        int rotNum = group->rotNum;
        float phi = bsList.at(i)->phi();
        float psi = bsList.at(i)->psi();
        char ss = bsList.at(i)->sscode;
        Phipsi pp(phi,psi);
        int ppIndex = ppLib.phipsiToIndex(&pp);
        char ppType = pp.regionAB();

        EnergyArray* eaS1 = new EnergyArray(rotNum);
        EnergyArray* eaRot = new EnergyArray(rotNum);
        EnergyArray* eaRef = new EnergyArray(rotNum);

        for(int j=0;j<rotNum;j++){
            Rotamer* rot = group->rotList.at(j);
            int aaType = rn.triToInt(rot->triName);
            float s1 = this->paList.at(i)->getScore(aaType);
            float ref = this->dp->ref[aaType];

            if(ss == 'H')
                ref = this->dp->refH[aaType];
            else if(ss == 'E')
                ref = this->dp->refE[aaType];
            else
                ref = this->dp->refC[aaType];

            float eRot;

            if(ppType == 'A')
                eRot = rotLibA->getNearestRotamerEnergy(rot, ppIndex);
            else
                eRot = rotLibB->getNearestRotamerEnergy(rot, ppIndex);

            eaS1->setEnergy(j, s1);
            eaRot->setEnergy(j, eRot);
            eaRef->setEnergy(j, ref);

        }
        this->s1EAList.push_back(eaS1);
        this->rotEAList.push_back(eaRot);
        this->refEAList.push_back(eaRef);
    }
}

void ScoringTemplate::loadPairwiseEnergy(){
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

        float sasaA = pair->siteA->data_[3];
        float sasaB = pair->siteB->data_[3];

        float meanSai = (sasaA+sasaB)/2;
        float vdwWT = dp->getVdwWeight(meanSai);


        ConformerGroup* groupA = cgList[idA];
        ConformerGroup* groupB = cgList[idB];

        EnergyMatrix* emS2 = new EnergyMatrix(groupA->rotNum, groupB->rotNum, idA, idB);
        EnergyMatrix* emVdw = new EnergyMatrix(groupA->rotNum, groupB->rotNum, idA, idB);

        for(int j=0;j<groupA->rotNum;j++){
            Conformer* confA = groupA->confList[j];
            int aaTypeA = rn.triToInt(confA->triName);
            for(int k=0;k<groupB->rotNum;k++){
                Conformer* confB = groupB->confList[k];
                int aaTypeB = rn.triToInt(confB->triName);
                float s2 = this->smList[i]->getValue(aaTypeA, aaTypeB);
                float vdw = pairEnergyScoring(confA, confB, seqSep, vdwWT, aec);
                float e = s2 + vdw;
                emS2->setEnergy(j,k, s2);
                emVdw->setEnergy(j,k, vdw);
            }
        }
        this->s2EMList.push_back(emS2);
        this->packEMList.push_back(emVdw);
    }

    for(int i=0;i<cgList.size();i++){
        delete cgList[i];
    }
}

void ScoringTemplate::updateResPairs() {
    /*
     * update bsPairs
     */
    this->bsPairs.clear();
    for(int i=0;i<this->resNum;i++){
        BackBoneSite* bsA = bsList[i];
        XYZ cbA = bsA->cbcrd();
        for(int j=i+1;j<this->resNum;j++){
            BackBoneSite* bsB = bsList[j];
            XYZ cbB = bsB->cbcrd();
            double d = cbA.distance(cbB);
            if(d < 8.0)
                bsPairs.push_back(new BackboneSitesPair(bsA, bsB));
            else if(d < 12 && contact(i,j,7.0))
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

void ScoringTemplate::loadS1S2(S1EnergyTable& s1Etable, S2MatrixFinder& s2Etable) {
    this->paList.clear();
    this->smList.clear();

    for(BackBoneSite* bs : bsList){
        AAProbabilityArray* pa = new AAProbabilityArray();
        s1Etable.getS1(*bs, pa);
        paList.push_back(pa);
    }

    for(BackboneSitesPair* pair : bsPairs){
        AAScoreMatrix* sm = new AAScoreMatrix();
        float d = pair->siteA->cbcrd().distance(pair->siteB->cbcrd());
        int sep = pair->seqSep;
        if(d < 8.0){
            s2Etable.getSM(pair, sm);
        }
        if(pair->siteA->sscode == 'E' && pair->siteB->sscode == 'E')
            sm->multiply(dp->wtS2EE);
        else if(pair->seqSep > 4)
            sm->multiply(dp->wtS2Nonlocal);
        else
            sm->multiply(dp->wtS2Local);


        smList.push_back(sm);
    }
}

void ScoringTemplate::setRotamerGroups() {
    int resNum = this->resList.size();
    vector<Rotamer*> rotList;
    for(int i=0;i<resNum;i++) {
        Rotamer* rot = resList[i]->natRotamer(&atLib);
        if(rot == NULL) {
            string s;
            rot = rotLibA->getAAGroup(this->resList[i]->triName)->rotList[0];
        }
        rotGroups[i]->clear();
        rotGroups[i]->addRotamer(rot);
    }
}

bool ScoringTemplate::contact(int posA, int posB, float cutoff){
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

float ScoringTemplate::clashScore(vector<Atom*>& listA, vector<Atom*>& listB, string typeB){
    /*
     * listA is backbone atoms
     * listB is sidechain atoms
     */
    float e = 0;
    XYZ zero;
    for(Atom* a : listA){
        string aName = a->getName();
        if(aName == "OXT")
            aName = "O";
        string uniqueNameA = "ALA-" + aName;
        //cout << "atomA: " << uniqueNameA << endl;
        AtomProperty* apA = atLib.getAtomProperty(uniqueNameA);
        for(Atom* b : listB){
            string uniqueNameB = typeB+"-"+b->getName();
            //cout << "atomB: " << uniqueNameB << endl;
            AtomProperty* apB = atLib.getAtomProperty(uniqueNameB);
            e += aec->getEnergy(a->distance(*b), apA->vdwRadius+apB->vdwRadius, -1.5, 1.0);
        }
    }
    return e;
}

void ScoringTemplate::deleteSc_BB_ClashedRotamers(float cutoff){
    for(int i=0;i<resNum;i++){
        //cout << "delete clash rotamer site: " << i+1 << endl;
        RotamerGroup* gp = rotGroups[i];
        BackBoneSite* bs = bsList[i];
        LocalFrame cs = getBackboneSiteLocalFrame(*bs);
        vector<int>* neighborList = neighborResidues[i];
        vector<Atom*> involvedAtoms;
        for(int j=0;j<neighborList->size();j++){
            RotamerGroup* gpNb = this->rotGroups[neighborList->at(j)];
            BackBoneSite* bsNb = this->bsList[neighborList->at(j)];
            int sep = abs(i - neighborList->at(j));
            if(sep <= 1)
                continue;
            involvedAtoms.push_back(new Atom("N", bsNb->ncrd()));
            involvedAtoms.push_back(new Atom("CA", bsNb->cacrd()));
            involvedAtoms.push_back(new Atom("C", bsNb->ccrd()));
            involvedAtoms.push_back(new Atom("O", bsNb->ocrd()));
        }

        float minClashScore = 9999;
        Rotamer* minClashRot = gp->rotList[0];

        for(int j=0;j<gp->rotNum;j++){
            //cout << "j=" << j << endl;
            Rotamer* rot = gp->rotList[j];
            //cout << "find rot: " << rot->rotName << endl;
            vector<Atom*> scAtoms;
            vector<XYZ> scTerns;

            rot->buildSidechain(cs, scTerns);
            //cout << "build sidechain" << endl;
            //cout << scTerns.size() << endl;
            for(int k=0;k<rot->atomNameList.size();k++){
                scAtoms.push_back(new Atom(rot->atomNameList[k], scTerns[k]));
            }
            //cout << "scAtoms: " << rot->triName << endl;
            float score = clashScore(involvedAtoms, scAtoms, rot->triName);


            for(int k=0;k<scAtoms.size();k++)
                delete scAtoms[k];

            if(score < minClashScore){
                minClashScore = score;
                minClashRot = rot;
            }

            if(score > cutoff){
                gp->deleteRotamer(j);
                j--;
            //    cout << "delete" << endl;
            }
        }

        if(gp->rotNum == 0){
            gp->addRotamer(minClashRot);
        }

        for(int k=0;k<involvedAtoms.size();k++)
            delete involvedAtoms[k];
    }
}



float pairEnergyScoring(Conformer* confA, Conformer* confB,int seqSep, float wt, AtomicEnergyCalcular* aec){

    float totalScSc = 0;
    float totalScBB = 0;

    int numScB = confB->scApList.size();
    int numScA = confA->scApList.size();

    /*
     * backboneA-sidechainB
    */
    string nameA, nameB;
    string triNameA, triNameB;
    XYZ *a, *b, *na, *nb;
    AtomProperty *apA, *apB;
    for(int i=0;i<4;i++){
        apA = confA->bbApList[i];
        a = confA->bbCoordList[i];
        for(int j=0;j<numScB;j++){
            b = confB->scCoordList[j];
            if(a->squaredDistance(*b) > 60.0) continue;
            nb = &(confB->normalVectorOfRing);
            apB = confB->scApList[j];
            float e = aec->getAtomEnergy(a, na, apA, b, nb, apB);
            if(e > 0)
                totalScBB += e;
            else
                totalScBB += e*wt;
        }
    }

    /*
     * sidechainA-backboneB
    */
    for(int i=0;i<numScA;i++){
        apA = confA->scApList[i];
        a = confA->scCoordList[i];
        na = &(confA->normalVectorOfRing);
        for(int j=0;j<4;j++){
            b = confB->bbCoordList[j];
            if(a->squaredDistance(*b) > 60.0) continue;
            apB = confB->bbApList[j];
            float e = aec->getAtomEnergy(a, na, apA, b, nb, apB);
            if(e > 0)
                totalScBB += e;
            else
                totalScBB += e*wt;
        }
    }

    /*
     * sidechainA-sidechainB
    */
    for(int i=0;i<numScA;i++){
        apA = confA->scApList[i];
        a = confA->scCoordList[i];
        na = &(confA->normalVectorOfRing);
        for(int j=0;j<numScB;j++){
            b = confB->scCoordList[j];
            if(a->squaredDistance(*b) > 60.0) continue;
            nb = &(confB->normalVectorOfRing);
            apB = confB->scApList[j];
            float e = aec->getAtomEnergy(a, na, apA, b, nb, apB);
            if(e > 0)
                totalScSc += e;
            else
                totalScSc += e*wt;
        }
    }

    if(seqSep == 1)
    {
        if(totalScSc > 0)
            return totalScSc;
        else
            return 0;
    }
    else
        return totalScBB + totalScSc;
}

ScoringTemplate::~ScoringTemplate() {
    delete rotLibA;
    delete rotLibB;
    delete aec;

    for(int i=0;i<s1EAList.size();i++){
        delete s1EAList.at(i);
    }

    for(int i=0;i<refEAList.size();i++){
        delete refEAList.at(i);
    }

    for(int i=0;i<rotEAList.size();i++){
        delete rotEAList.at(i);
    }

    for(int i=0;i<rotGroups.size();i++){
        delete rotGroups.at(i);
    }

    for(int i=0;i<s2EMList.size();i++){
        delete s2EMList.at(i);
    }

    for(int i=0;i<packEMList.size();i++){
        delete packEMList.at(i);
    }

    for(int i=0;i<paList.size();i++){
        delete paList.at(i);
    }

    for(int i=0;i<smList.size();i++){
        delete smList.at(i);
    }

    for(int i=0;i<involvedPairs.size();i++){
        delete involvedPairs.at(i);
    }

    for(int i=0;i<neighborResidues.size();i++){
        delete neighborResidues.at(i);
    }
}

} /* namespace NSPtest */

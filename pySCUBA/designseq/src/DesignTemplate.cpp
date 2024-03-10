/*
 * DesignTemplate.cpp
 *
 *  Created on: 2017��12��2��
 *      Author: notxp
 */

#include "designseq/DesignTemplate.h"

namespace NSPdesignseq {

DesignTemplate::DesignTemplate(vector<BackBoneSite*>& bsList, DesignParameters* dp, 
    S1EnergyTable& s1Etable, S2MatrixFinder& s2Etable, double vdwres) {

    backboneSiteListUpdateSasa(bsList);
    this->s1ET = &s1Etable;
    this->resNum = bsList.size();
    this->aec = new AtomicEnergyCalcular(dp, vdwres);
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

    SeqProfile prof;


    getSeqProfile(s1Etable, s2Etable, &prof);

    addProfileRotamerGroups(&prof, 0.3);
    updateResPairs();

    deleteSc_BB_ClashedRotamers(6.0);
    updateResPairs();
    this->natRefEnergy = 0.0;
}

DesignTemplate::DesignTemplate(vector<BackBoneSite*>& bsList, string& resFile, DesignParameters* dp,
    S1EnergyTable& s1Etable, S2MatrixFinder& s2Etable, double vdwres) {
    backboneSiteListUpdateSasa(bsList);
    this->s1ET = &s1Etable;
    this->resNum = bsList.size();
    this->aec = new AtomicEnergyCalcular(dp,vdwres);
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

    SeqProfile prof;
    getSeqProfile(s1Etable, s2Etable, &prof);
    addRotamerGroupsFromResFile(resFile, &prof, 0.3);
    updateResPairs();
    deleteSc_BB_ClashedRotamers(6.0);
    updateResPairs();
    this->natRefEnergy = 0.0;
}

void DesignTemplate::loadSingleResidueEnergy(){

    PhipsiLib ppLib;
    ecList.clear();
    for(int i=0;i<resNum;i++){
        ecList.push_back(EnergyComponents());
        RotamerGroup* group = rotGroups.at(i);
        int rotNum = group->rotNum;
        cout << "res: " << i+1 << " rotNum: " << rotNum << endl;
        float phi = bsList.at(i)->phi();
        float psi = bsList.at(i)->psi();
        char ss = bsList.at(i)->sscode;
        double sai = bsList.at(i)->data_[3];
        Phipsi pp(phi,psi);
        int ppIndex = ppLib.phipsiToIndex(&pp);
        char ppType = pp.regionAB();
        EnergyArray* ea = new EnergyArray(rotNum);
        for(int j=0;j<rotNum;j++){
        //    cout << "rotID: " << j << endl;
            Rotamer* rot = group->rotList.at(j);
        //    cout << "rotType: " << rot->rotName     << endl;
            int aaType = rn.triToInt(rot->triName);
        //    cout << "aaType: " << aaType << endl;
            float s1 = this->paList.at(i)->getScore(aaType);
        //    cout << "s1: " << s1 << endl;
/*            float ref = this->dp->ref[aaType];

            if(ss == 'H')
                ref = this->dp->refH[aaType];
            else if(ss == 'E')
                ref = this->dp->refE[aaType];
            else
                ref = this->dp->refC[aaType];
*/
            float ref = this->dp->getRefELocal(ss, sai, rot->triName);
        //    cout << "ref: " << ref << endl;
            float eRot;

        //    float eRotInterpolation;
            if(ppType == 'A')
                eRot = rotLibA->getRotamerEnergy(rot->rotName, ppIndex);
            else
                eRot = rotLibB->getRotamerEnergy(rot->rotName, ppIndex);

            /*
            if(ppType == 'A')
                eRotInterpolation = rotLibA->getRotamerEnergy(rot->rotName, &pp);
            else
                eRotInterpolation = rotLibB->getRotamerEnergy(rot->rotName, &pp);

            printf("eRot Interpolation: %6.4f %6.4f\n", eRot, eRotInterpolation );
             */
        //    cout << "rotEnergy: " << eRot << endl;
            ea->setEnergy(j, s1+ref+eRot);
            ecList.back().s1.push_back(s1);
            ecList.back().ref.push_back(ref);
            ecList.back().rot.push_back(eRot);
        }
        this->eaList.push_back(ea);
    }
}



void DesignTemplate::loadPairwiseEnergy(){

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
    //    cout << "groupA" << endl;
        ConformerGroup* groupB = cgList[idB];
    //    cout << "groupB" << endl;


//        RotamerGroup* groupA = rotGroups.at(idA);
//        RotamerGroup* groupB = rotGroups.at(idB);
        EnergyMatrix* em = new EnergyMatrix(groupA->rotNum, groupB->rotNum, idA, idB);

        for(int j=0;j<groupA->rotNum;j++){
            Conformer* confA = groupA->confList[j];
        //    cout << "confA " << j << endl;
            int aaTypeA = rn.triToInt(confA->triName);
            for(int k=0;k<groupB->rotNum;k++){
                Conformer* confB = groupB->confList[k];
        //        cout << "confB " << k <<  endl;
                int aaTypeB = rn.triToInt(confB->triName);
                float s2 = this->smList[i]->getValue(aaTypeA, aaTypeB);

                float vdw = pairEnergy(confA, confB, seqSep, vdwWT, aec);

                float e = s2 + vdw;
                if(e > 10000 || e < -1000)
                {
                    cerr << "aaType: " << aaTypeA << " " << aaTypeB << endl;
                    cerr << "s2: " << s2 << endl;
                    cerr << "vdw: " << vdw << endl;
                    cerr << "calculate energy error: pos-" << idA << " pos-" << idB << " rotA-" << j << " rotB-" << k << endl;
                    exit(1);
                }
                em->setEnergy(j,k, s2 + vdw);
            }
        }

        for(int j=0;j<groupA->rotNum;j++){
            for(int k=0;k<groupB->rotNum;k++){
                float e = em->getEnergy(j,k);
                if(e > 10000 || e < -1000)
                {
                    cerr << "assign energy error: pos-" << idA << " pos-" << idB << " rotA-" << j << " rotB-" << k << endl;
                    exit(1);
                }
            }
        }

        this->emList.push_back(em);
    }

    for(int i=0;i<cgList.size();i++){
        delete cgList[i];
    }
}

void DesignTemplate::updateResPairs(){
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
            else if(d < 12 && contact(i,j,6.0))
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

void DesignTemplate::loadS1S2(S1EnergyTable& s1Etable, S2MatrixFinder& s2Etable){
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

void DesignTemplate::updateNativeSequenceEnergy(double eNat){
    //ResName rn;
    for(int i=0;i<resNum;i++){
        RotamerGroup* group = rotGroups[i];
        int rotNum = group->rotNum;
        EnergyArray* ea = eaList[i];
        for(int j=0;j<rotNum;j++){
            Rotamer* rot = group->rotList.at(j);
            int aaType = rn.triToInt(rot->triName);
            int natAAType = rn.triToInt(bsList[i]->resname);
            if(aaType == natAAType)
                ea->setEnergy(j, ea->getEnergy(j)+eNat);
        }
    }
}

void DesignTemplate::getSeqProfile(S1EnergyTable& s1Etable, S2MatrixFinder& s2Etable, SeqProfile* result){

    AADesignTemplate aadt(this->bsList, &s1Etable, &s2Etable);
    aadt.getProfile(result);
}

bool DesignTemplate::contact(int posA, int posB, float cutoff){
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

void DesignTemplate::setRotamerGroups(vector<RotamerGroup*>& group){
    for(int pos=0;pos<this->resNum;pos++){
        rotGroups[pos]->clear();
        int rotNum = group[pos]->rotNum;
        for(int j=0;j<rotNum;j++){
            rotGroups[pos]->addRotamer(group[pos]->rotList[j]);
        }
    }
}

void DesignTemplate::addProfileRotamerGroups(SeqProfile* prof, float cutoff){
/*
 * suggest cutoff 0.2
 * Cysteine is not allowed
 */

    //cout << "print profile: " << endl;
    //prof->printProfile();

    PhipsiLib ppLib;

    for(int pos=0;pos<this->resNum;pos++){
        float phi = bsList.at(pos)->phi();
        float psi = bsList.at(pos)->psi();
        AAProbabilityArray pa;
        this->s1ET->getS1(*bsList[pos], &pa);

        Phipsi pp(phi, psi);
        int ppType = ppLib.phipsiToIndex(&pp);
        RotamerLib* rotLib;
        if(pp.regionAB() == 'A')
            rotLib = rotLibA;
        else
            rotLib = rotLibB;
        RotamerGroup* gp = rotGroups[pos];
        for(int i=0;i<20;i++){
            if(i==1) continue;
            string triName = rn.intToTri(i);

            float s1relP = pa.getRelP(i);
            float relP = prof->getRelP(pos, i);
            if((relP > cutoff && s1relP > cutoff) || i == 5){
                RotamerGroup* gp0 = rotLib->getAAGroup(triName);
                for(Rotamer* rot : gp0->rotList){
                    float eRot = rotLib->getRotamerEnergy(rot->rotName, ppType);
                    if(eRot < 4)
                        gp->addRotamer(rot);
                }
            }
        }
    //    cout << "5: " << endl;
    }
}

void DesignTemplate::addRotamerGroupsFromResFile(string& resFile, SeqProfile* prof, float cutoff){
    ifstream rf;
    rf.open(resFile.c_str(),ios::in);
    if(!rf.is_open()){
        cout << "fail to open file " << resFile << endl;
        exit(1);
    }
    string s;

    cout << "selected amino acid choice from resFile" << endl;

    vector<string> posAAChoice;
    for(int i=0;i<resNum;i++){
        posAAChoice.push_back("");
    }

    string allButCys = "ADEFGHIKLMNPQRSTVWY";
    string all = "ACDEFGHIKLMNPQRSTVWY";
    string defaultValue = "nat";

    vector<string> spt;
    while(getline(rf,s)){
        if (s.empty()) {
            continue;
        }
        if (s[0] == '#') {
            continue;
        }
        splitString(s," ",&spt);
        if(spt.size() == 2 && spt[0] == "default") {
            if(spt[1] == "nat") {
                for(int i=0;i<resNum;i++){
                    string triName = this->bsList[i]->resname;
                    char xx[2];
                    xx[0] = rn.triToSin(triName);
                    xx[1] = '\0';
                    posAAChoice[i] = string(xx);
                }
            }
            else if(spt[1] == "all") {
                for(int i=0;i<resNum;i++){
                    posAAChoice[i] = all;
                }
            }
            else if(spt[1] == "allButCys") {
                for(int i=0;i<resNum;i++){
                    posAAChoice[i] = allButCys;
                }
            }
            else {
                for(int i=0;i<resNum;i++){
                    posAAChoice[i] = spt[1];
                }
            }
            continue;
        }
        if(spt.size() != 3){
            cout << "invalid res choice: " << s << endl;
            exit(0);
        }

        string key = spt[0] + spt[1];
        if(chainIDResIDToIndex.find(key) == chainIDResIDToIndex.end()){
            cout << "invalid resID: " << s << endl;
            exit(0);
        }

        int index = chainIDResIDToIndex[key];
        if(spt[2] == "all")
            posAAChoice[index] = all;
        else if(spt[2] == "allButCys")
            posAAChoice[index] = allButCys;
        else if(spt[2] == "nat"){
            string triName = this->bsList.at(index)->resname;
            char xx[2];
            xx[0] = rn.triToSin(triName);
            xx[1] = '\0';
            posAAChoice[index] = string(xx);
        }
        else
            posAAChoice[index] = spt[2];
    }

    rf.close();

    PhipsiLib ppLib;

    for(int pos=0;pos<this->resNum;pos++){
        float phi = bsList.at(pos)->phi();
        float psi = bsList.at(pos)->psi();
        AAProbabilityArray pa;
        this->s1ET->getS1(*bsList[pos], &pa);

        Phipsi pp(phi, psi);
        int ppType = ppLib.phipsiToIndex(&pp);
        RotamerLib* rotLib;
        if(pp.regionAB() == 'A')
            rotLib = rotLibA;
        else
            rotLib = rotLibB;
        RotamerGroup* gp = rotGroups[pos];
        string aaChoice = posAAChoice[pos];
        if(aaChoice == ""){
            for(int i=0;i<20;i++){
                if(i==1) continue;
                string triName = rn.intToTri(i);

                float s1relP = pa.getRelP(i);
                float relP = prof->getRelP(pos, i);
                if((relP > cutoff && s1relP > cutoff) || i == 5){
                    RotamerGroup* gp0 = rotLib->getAAGroup(triName);
                    for(Rotamer* rot : gp0->rotList){
                        float eRot = rotLib->getRotamerEnergy(rot->rotName, ppType);
                        if(eRot < 4)
                            gp->addRotamer(rot);
                    }
                }
            }
        }
        else{
            cout << "pos: " << pos << " choice:" << aaChoice << endl;
            std::vector<float> relP_n(20,0.0);
            std::vector<float> s1relP_n(20,0.0);
            float relP_t = 0.0;
            float s1relP_t = 0.0;
            float r0 = 0.0;
            float s1r0 = 0.0;
            for(int i = 0; i < 20; i++){
                r0 += prof->getRelP(pos, i) + 1.e-8;
                s1r0 += pa.getRelP(i) + 1.e-8;
            }
            for(int i = 0; i < aaChoice.length(); i++){
                int aaTypeInt = rn.sinToInt(aaChoice[i]);
                s1relP_t += pa.getRelP(aaTypeInt) + 1.e-8;
                relP_t += prof->getRelP(pos, aaTypeInt) + 1.e-8;
            }
            for(int i = 0; i < aaChoice.length(); i++){
                 int aaTypeInt = rn.sinToInt(aaChoice[i]);
                 s1relP_n[aaTypeInt] = (pa.getRelP(aaTypeInt) + 1.e-8) * s1r0 / s1relP_t;
                 relP_n[aaTypeInt] = (prof->getRelP(pos, aaTypeInt) + 1.e-8) * r0 / relP_t;
             }
            for(int i=0;i<aaChoice.length();i++){
                int aaTypeInt = rn.sinToInt(aaChoice[i]);
                string aaTypeTri = rn.sinToTri(aaChoice[i]);
                if(aaTypeInt < 0 || aaTypeInt > 19){
                    cout << "invalid aa choice: " << aaChoice << endl;
                    exit(0);
                }
                RotamerGroup* gp0 = rotLib->getAAGroup(aaTypeTri);
                if((relP_n[aaTypeInt] <= cutoff || s1relP_n[aaTypeInt] <=cutoff)
                        && aaTypeInt != 5) continue;
                for(Rotamer* rot : gp0->rotList){
                    float eRot = rotLib->getRotamerEnergy(rot->rotName, ppType);
                    if(eRot < 4)
                        gp->addRotamer(rot);
                }
            }
        }
    }
}


float DesignTemplate::clashScore(vector<Atom*>& listA, vector<Atom*>& listB, string typeB){
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

void DesignTemplate::deleteSc_BB_ClashedRotamers(float cutoff){
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

float pairEnergy(Conformer* confA, Conformer* confB,int seqSep, float wt, AtomicEnergyCalcular* aec){

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

float pairEnergyHB(Conformer* confA, Conformer* confB,int seqSep, float wt, AtomicEnergyCalcular* aec){

    float totalScSc = 0;
    float totalScBB = 0;

    /*
     * backboneA-sidechainB
    */
    string nameA, nameB;
    string triNameA, triNameB;
    XYZ *a, *b, *na, *nb;
    AtomProperty *apA, *apB;

    vector<PolarAtom*> paList = confA->bbPolarList;
    vector<PolarAtom*> pbList = confB->scPolarList;
    for(int i=0;i<paList.size();i++){
        PolarAtom* pa = paList[i];
        for(int j=0;j<pbList.size();j++){
            PolarAtom* pb = pbList[j];
            double e = aec->hbEnergy(pa, pb);
            totalScBB += e;
        }
    }

    /*
     * sidechainA-backboneB
    */
    paList = confA->scPolarList;
    pbList = confB->bbPolarList;
    for(int i=0;i<paList.size();i++){
        PolarAtom* pa = paList[i];
        for(int j=0;j<pbList.size();j++){
            PolarAtom* pb = pbList[j];
            double e = aec->hbEnergy(pa, pb);
            totalScBB += e;
        }
    }

    /*
     * sidechainA-sidechainB
    */
    paList = confA->scPolarList;
    pbList = confB->scPolarList;
    for(int i=0;i<paList.size();i++){
        PolarAtom* pa = paList[i];
        for(int j=0;j<pbList.size();j++){
            PolarAtom* pb = pbList[j];
            double e = aec->hbEnergy(pa, pb);
            totalScSc += e;
        }
    }

    if(seqSep == 1)
    {
        return 0;
    }
    else
        return totalScBB + totalScSc;
}

void DesignTemplate::printRotChoice(){
    for(int i=0;i<this->resNum;i++){
        int rotNum = this->rotGroups.at(i)->rotNum;
        printf("%-4d %3d",i+1, rotNum);
        for(int j=0;j<rotNum;j++){
            printf(" %-6s", this->rotGroups.at(i)->rotList.at(j)->rotName.c_str());
        }
        cout << endl;
    }
}

DesignTemplate::~DesignTemplate() {

    delete rotLibA;
    delete rotLibB;
    delete aec;

    for(int i=0;i<eaList.size();i++){
        delete eaList.at(i);
    }

    for(int i=0;i<rotGroups.size();i++){
        delete rotGroups.at(i);
    }

    for(int i=0;i<emList.size();i++){
        delete emList.at(i);
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
    // TODO Auto-generated destructor stub
}



} /* namespace NSPdesignseq */

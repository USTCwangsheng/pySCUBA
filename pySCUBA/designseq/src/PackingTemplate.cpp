/*
 * PackingTemplate.cpp
 *
 *  Created on: 2018Äê4ÔÂ16ÈÕ
 *      Author: notxp
 */

#include "designseq/PackingTemplate.h"

namespace NSPdesignseq {

PackingTemplate::PackingTemplate(PDB* pdb, DesignParameters* para) {

    this->pdb = pdb;
    float rotECutoff = 5.0;
    this->para = para;

    this->rotLibA = new RotamerLib("A1");
    this->rotLibB = new RotamerLib("B1");
    this->atLib = new AtomLib();
    this->ppLib = new PhipsiLib();
    this->ec = new AtomicEnergyCalcular(para);

    ResName rn;
    size_t i,j;
    vector<Residue*> pdbList0 = pdb->getResList();
    vector<Residue*> pdbResList;
    for(i=0;i<pdbList0.size();i++) {
        Residue* res = pdbList0[i];
        if(!res->hasThreeCoreAtoms() || !res->hasAtom("O"))
        {
            continue;
        }
        pdbResList.push_back(res);
    }


    StructureInfo si(pdb);
    si.updateTorsion();
    SasaPSD rsp;
    si.updateSAI(&rsp);

    int x = si.getResNum();

    int posIndex = 0;
    for(i=0;i<pdbResList.size();i++)
    {
        Residue* res = pdbResList.at(i);

        if(!res->hasThreeCoreAtoms() || !res->hasAtom("O"))
        {
            cout << "lack core atom: " << res->resID << endl;
            continue;
        }
        if(res->triName == "MSE")
            res->triName = "MET";
        if(!rn.isStandardAminoAcid(res->triName))
        {
            cout << "invalid residue: " << res->resID << " " << res->triName << endl;
            continue;
        }

        BackBoneSite* bs = new BackBoneSite();

    //    cout << "resID: " << res->resID    << endl;
        bs->pdbid = pdb->getPDBID();
        bs->chainid = res->chainID;
        bs->strresid = si.getResID(i);
        bs->resname = res->getType();
        bs->resseq = i;

        bs->sscode = si.getSS(i);
        bs->data_[0] = si.getPhi(i);
        bs->data_[1] = si.getPsi(i);
        bs->data_[2] = si.getOmg(i);
        bs->data_[3] = si.getSai(i);

    //    cout << "STRUCTUREINFO: PC2BBSITES: " << bs.resid << " " << bs.resname << " " << bs.data_[0] << " " << bs.data_[1] << " " << bs.data_[3] << endl;
        Atom* aN = res->getAtom("N");
        Atom* aCA = res->getAtom("CA");
        Atom* aC = res->getAtom("C");
        Atom* aO = res->getAtom("O");
        XYZ N = aN->getCoord();
        XYZ CA = aCA->getCoord();
        XYZ C = aC->getCoord();
        XYZ O = aO->getCoord();
        bs->data_[4] = N[0];
        bs->data_[5] = N[1];
        bs->data_[6] = N[2];
        bs->data_[7] = CA[0];
        bs->data_[8] = CA[1];
        bs->data_[9] = CA[2];
        bs->data_[10] = C[0];
        bs->data_[11] = C[1];
        bs->data_[12] = C[2];
        bs->data_[13] = O[0];
        bs->data_[14] = O[1];
        bs->data_[15] = O[2];



        bsList.push_back(bs);
        RotamerGroup* gp = new RotamerGroup();
        rotGroups.push_back(gp);
    }

    this->resNum = bsList.size();

    updateRotamerGroups();
    updateResPairs();
    updateEAEM();
}

void PackingTemplate::updateRotamerGroups(){

    for(int i=0;i<this->resNum;i++){

        string triName = this->bsList[i]->resname;
        if(triName != "CYS") continue;

        XYZ cbA = bsList[i]->cbcrd();

        float phi = bsList.at(i)->phi();
        float psi = bsList.at(i)->psi();

        Phipsi pp(phi, psi);
        int ppType = ppLib->phipsiToIndex(&pp);

        RotamerLib* rotLib;
        if(pp.regionAB() == 'A')
            rotLib = rotLibA;
        else
            rotLib = rotLibB;

        RotamerGroup* gpI = rotLib->getAAGroup(triName);
        ConformerGroup* cgI = new ConformerGroup(gpI, bsList[i], atLib);

        for(int j=i+1;j<this->resNum;j++){

            triName = this->bsList[j]->resname;
            if(triName != "CYS") continue;
            XYZ cbB = bsList[i]->cbcrd();
            if(cbA.distance(cbB) > 6) continue;
            phi = bsList.at(j)->phi();
            psi = bsList.at(j)->psi();

            Phipsi pp2(phi, psi);
            int ppType2 = ppLib->phipsiToIndex(&pp2);

            if(pp2.regionAB() == 'A')
                rotLib = rotLibA;
            else
                rotLib = rotLibB;

            RotamerGroup* gpJ = rotLib->getAAGroup(triName);
            ConformerGroup* cgJ = new ConformerGroup(gpJ, bsList[j], atLib);

            double minE = 999;
            int minEIndexI = -1;
            int minEIndexJ = -1;
            for(int ci=0;ci<cgI->rotNum;ci++){
                Conformer* confA = cgI->confList[ci];
                for(int cj=0;cj<cgJ->rotNum;cj++){
                    Conformer* confB = cgJ->confList[cj];
                    XYZ* cbI = confA->scCoordList[0];
                    XYZ* sgI = confA->scCoordList[1];
                    XYZ* cbJ = confB->scCoordList[0];
                    XYZ* sgJ = confB->scCoordList[1];
                    float dSS = sgI->distance(*sgJ);
                    float ang1 = angleX(*cbI, *sgI, *sgJ);
                    float ang2 = angleX(*sgI, *sgJ, *cbJ);
                    float dihed = dihedral(*cbI, *sgI, *sgJ, *cbJ);

                    float eSS = (dSS-2.055)*(dSS-2.055);
                    eSS += 0.00063*(ang1-105.2)*(ang1-105.2);
                    eSS += 0.00063*(ang2-105.2)*(ang2-105.2);
                    eSS += 0.00018*(abs(dihed)-89.7)*(abs(dihed)-89.7);
                    if(eSS < minE){
                        minE = eSS;
                        minEIndexI = ci;
                        minEIndexJ = cj;
                    }
                }
            }
            if(minE < 0.3){
                rotGroups[i]->addRotamer(gpI->rotList[minEIndexI]);
                rotGroups[j]->addRotamer(gpJ->rotList[minEIndexJ]);
            }
        }
    }


    for(int pos=0;pos<this->resNum;pos++){

        float phi = bsList.at(pos)->phi();
        float psi = bsList.at(pos)->psi();

        Phipsi pp(phi, psi);
        int ppType = ppLib->phipsiToIndex(&pp);

        RotamerLib* rotLib;
        if(pp.regionAB() == 'A')
            rotLib = rotLibA;
        else
            rotLib = rotLibB;
        RotamerGroup* gp = rotGroups[pos];
        if(gp->rotNum == 1)
            continue;
        string triName = bsList[pos]->resname;
        if(triName == "MSE")
            triName = "MET";

        RotamerGroup* gp0 = rotLib->getAAGroup(triName);
        for(Rotamer* rot : gp0->rotList){
            float eRot = rotLib->getRotamerEnergy(rot->rotName, ppType);
            if(eRot < 4.5)
                gp->addRotamer(rot);
        }
    }
}

bool PackingTemplate::contact(int posA, int posB, float cutoff){
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


void PackingTemplate::updateResPairs(){
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
            {
                bsPairs.push_back(new BackboneSitesPair(bsA, bsB));

            }
            else if(d < 12 && contact(i,j,6.0))
                bsPairs.push_back(new BackboneSitesPair(bsA, bsB));
        }
    }


    /*
     * update involvedPairs and neighbor residues
     */


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



void PackingTemplate::updateEAEM()
{
    PhipsiLib ppLib;
    for(int i=0;i<resNum;i++){
        RotamerGroup* group = rotGroups[i];
        int rotNum = group->rotNum;
    //    cout << "res: " << i+1 << " rotNum: " << rotNum << endl;
        float phi = bsList[i]->phi();
        float psi = bsList[i]->psi();
        Phipsi pp(phi,psi);
        int ppIndex = ppLib.phipsiToIndex(&pp);
        char ppType = pp.regionAB();
        EnergyArray* ea = new EnergyArray(rotNum);
        for(int j=0;j<rotNum;j++){
            Rotamer* rot = group->rotList[j];
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

    vector<ConformerGroup*> cgList;
    cgList.push_back(new ConformerGroup(rotGroups[0], bsList[0], atLib));

    for(int i=1;i<resNum;i++){
        cgList.push_back(new ConformerGroup(rotGroups[i], bsList[i], bsList[i-1], atLib));
    }


    for(int i=0;i<this->bsPairs.size();i++){

        BackboneSitesPair* pair = this->bsPairs[i];
        int idA = pair->siteA->resseq;
        int idB = pair->siteB->resseq;
        int seqSep = abs(idA-idB);
        if(seqSep > 5) seqSep = 5;
        if(pair->siteA->chainid != pair->siteB->chainid)
            seqSep = 5;

        ConformerGroup* groupA = cgList[idA];
        ConformerGroup* groupB = cgList[idB];


        EnergyMatrix* em = new EnergyMatrix(groupA->rotNum, groupB->rotNum, idA, idB);

        for(int j=0;j<groupA->rotNum;j++){
            Conformer* confA = groupA->confList[j];
            int aaTypeA = rn.triToInt(confA->triName);
            for(int k=0;k<groupB->rotNum;k++){
                Conformer* confB = groupB->confList[k];
                int aaTypeB = rn.triToInt(confB->triName);

                float sasaA = pair->siteA->data_[3];
                float sasaB = pair->siteB->data_[3];

                float meanSai = (sasaA+sasaB)/2;
                float vdwWT = this->para->getVdwWeight(meanSai);
                float vdw = pairEnergyPack(confA, confB, seqSep, vdwWT, ec);
                //float eHb = pairHBEnergyPack(confA, confB, ec);

                //float e = vdw;

                em->setEnergy(j,k, vdw);
            }
        }

        this->emList.push_back(em);
    }

    for(int i=0;i<cgList.size();i++){
        delete cgList[i];
    }

}

void PackingTemplate::printDetail() {
    int posNum = this->bsList.size();
    int pairNum = this->bsPairs.size();

    cout << "Single site Energy: " << endl;
    for(int i=0;i<posNum;i++) {
        RotamerGroup* gp = this->rotGroups.at(i);
        int choiceNum = gp->rotNum;
        for(int j=0;j<choiceNum;j++) {
            Rotamer* rot = gp->rotList[j];
            double ene = this->eaList[i]->getEnergy(j);
            printf("res: %-2d %-7s %8.3f\n",i, rot->rotName.c_str(), ene);
        }
    }

    cout << "Pairwise Energy" << endl;
    for(int i=0;i<pairNum;i++) {
        BackboneSitesPair* pair = bsPairs[i];
        int numA = this->emList[i]->choiceANum;
        int numB = this->emList[i]->choiceBNum;
        int posA = this->emList[i]->posA;
        int posB = this->emList[i]->posB;
        for(int j=0;j<numA;j++) {
            Rotamer* rotA = this->rotGroups[posA]->rotList[j];
            for(int k=0;k<numB;k++) {
                Rotamer* rotB = this->rotGroups[posB]->rotList[k];
                double ene = this->emList[i]->getEnergy(j,k);
                printf("pair: %-2d %-2d %-7s %-7s %8.3f\n",this->emList[i]->posA, this->emList[i]->posB, rotA->rotName.c_str(), rotB->rotName.c_str(), ene);

            }
        }
    }
}

void PackingTemplate::testUbq() {
    BackBoneSite* bsA = this->bsList[14];
    Rotamer* rotA = this->rotGroups[14]->rotList[0];

    BackBoneSite* bsB = this->bsList[25];
    Rotamer* rotB = this->rotGroups[25]->rotList[0];

    Conformer* confA = new Conformer(rotA, bsA, this->atLib);
    Conformer* confB = new Conformer(rotB, bsB, this->atLib);



}

float pairEnergyPack(Conformer* confA, Conformer* confB,int seqSep, float wt, AtomicEnergyCalcular* aec){


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
            float e = aec->getAtomEnergyPacking(a, na, apA, b, nb, apB);

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
            float e = aec->getAtomEnergyPacking(a, na, apA, b, nb, apB);

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

            float e = aec->getAtomEnergyPacking(a, na, apA, b, nb, apB);

            if(e > 0)
                totalScSc += e;
            else
                totalScSc += e*wt;
        }
    }

    float eDS = 0;
    if(confA->triName == "CYS" && confB->triName == "CYS"){
        XYZ* cbI = confA->scCoordList[0];
        XYZ* sgI = confA->scCoordList[1];
        XYZ* cbJ = confB->scCoordList[0];
        XYZ* sgJ = confB->scCoordList[1];
        float dSS = sgI->distance(*sgJ);
        float ang1 = angleX(*cbI, *sgI, *sgJ);
        float ang2 = angleX(*sgI, *sgJ, *cbJ);
        float dihed = dihedral(*cbI, *sgI, *sgJ, *cbJ);

        float eSS = (dSS-2.055)*(dSS-2.055);
        eSS += 0.00063*(ang1-105.2)*(ang1-105.2);
        eSS += 0.00063*(ang2-105.2)*(ang2-105.2);
        eSS += 0.00018*(abs(dihed)-89.7)*(abs(dihed)-89.7);
        eDS = (eSS - 0.5)*20.0;
        //cout << "dsWT " << aec->dsWT << " " << eDS<< endl;
    }

    if(eDS != 0 && eDS < totalScSc)
        totalScSc = eDS;

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


float pairHBEnergyPack(Conformer* confA, Conformer* confB, int index, AtomicEnergyCalcular* aec){

    float hbEnergy = 0;
    int numScB = confB->scApList.size();
    int numScA = confA->scApList.size();

    float eBbSc = 0;
    float eScBb = 0;
    float eScSc = 0;

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
            if(a->squaredDistance(*b) > 49.0) continue;
            nb = &(confB->normalVectorOfRing);
            apB = confB->scApList[j];
            float e = aec->getHBAtomEnergyPacking(a, na, apA, b, nb, apB);
            eBbSc += e;
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
            if(a->squaredDistance(*b) > 49.0) continue;
            apB = confB->bbApList[j];
            float e = aec->getHBAtomEnergyPacking(a, na, apA, b, nb, apB);
            eScBb += e;
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
            if(a->squaredDistance(*b) > 49.0) continue;
            nb = &(confB->normalVectorOfRing);
            apB = confB->scApList[j];

            float e = aec->getHBAtomEnergyPacking(a, na, apA, b, nb, apB);

            if(apA->isPolar && apB->isPolar) {
                float d = a->distance(*b);
            }
            eScSc += e;
        }
    }


    if(index == 0)
        return eScBb + eScSc;
    else
        return eBbSc + eScSc;
}



PackingTemplate::~PackingTemplate() {
    // TODO Auto-generated destructor stub
    size_t i;
    delete rotLibA;
    delete rotLibB;
    delete atLib;
    delete ppLib;
    delete ec;

    for(int i=0;i<eaList.size();i++){
        delete eaList.at(i);
    }

    for(int i=0;i<rotGroups.size();i++){
        delete rotGroups.at(i);
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


}

} /* namespace packing */

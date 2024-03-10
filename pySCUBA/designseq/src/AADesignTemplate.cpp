/*
 * AADesignTemplate.cpp
 *
 *  Created on: 2017Äê12ÔÂ15ÈÕ
 *      Author: notxp
 */

#include "designseq/AADesignTemplate.h"

namespace NSPdesignseq {

AADesignTemplate::AADesignTemplate(vector<BackBoneSite*>& bsList) {
    for(BackBoneSite* bs : bsList){
        this->bsList.push_back(bs);

    }

    XYZ localCB(-0.947 , 0.001  ,1.197);
    int n = bsList.size();
    for(int i=0;i<n;i++){
        BackBoneSite* bsA = bsList[i];
        LocalFrame csA = getBackboneSiteLocalFrame(*bsA);
        XYZ cbA = csA.local2globalcrd(localCB);
        for(int j=i+1;j<n;j++){
            BackBoneSite* bsB = bsList[j];
            LocalFrame csB = getBackboneSiteLocalFrame(*bsB);
            XYZ cbB = csB.local2globalcrd(localCB);
            if(cbA.distance(cbB) < 8.0)
            {
                bsPairs.push_back(new BackboneSitesPair(bsA, bsB));
            }
        }
    }
}

AADesignTemplate::AADesignTemplate(vector<BackBoneSite*>& bsList, S1EnergyTable* s1Etable, S2MatrixFinder* s2Etable) {

    for(BackBoneSite* bs : bsList){
        this->bsList.push_back(bs);

    }

    int n = bsList.size();
    for(int i=0;i<n;i++){
        BackBoneSite* bsA = bsList[i];
        XYZ cbA = bsA->cbcrd();
        for(int j=i+1;j<n;j++){
            BackBoneSite* bsB = bsList[j];
            XYZ cbB = bsB->cbcrd();
            if(cbA.distance(cbB) < 8.0)
            {
                bsPairs.push_back(new BackboneSitesPair(bsA, bsB));
            }
        }
    }
    loadS1S2(s1Etable, s2Etable);
}



void AADesignTemplate::loadS1S2(string& paFile, string& smFile) {
    ifstream input;
    input.open(paFile.c_str(), ios::in);
    //cout << "open paFile: " << paFile << endl;
    if(! input.is_open()){
        cout << "fail to open file: " << paFile << endl;
    }
    string s;
    vector<string> spt;
    while(getline(input,s)){
        splitString(s, " ", &spt);
        double p[20];
        for(int i=0;i<20;i++){
            p[i] = exp(-1.0*atof(spt[i].c_str()));
        }
        AAProbabilityArray* pa = new AAProbabilityArray(p);
        this->paList.push_back(pa);
    }

    input.close();


    input.open(smFile.c_str(), ios::in);
    //cout << "open smFile: " << smFile << endl;
    if(!input.is_open()){
        cout << "fail to open file: " << smFile << endl;
    }

    float matrix[20][20];

    int lineID = -2;
    while(getline(input, s)){
        if(lineID < 0) {
            lineID++;
            continue;
        }
        splitString(s, " ", &spt);
        for(int i=0;i<20;i++){
            matrix[lineID][i] = atof(spt[i].c_str());
        }
        lineID++;
        if(lineID == 20){
            AAScoreMatrix* sm = new AAScoreMatrix();
            sm->initValue(matrix);
            int seqSep = bsPairs[this->smList.size()]->seqSep;
            double wt0 = 0.8;
        //    double wt = 0.5;
            if(seqSep < 5)
                wt0 = 0.5;


            sm->multiply(wt0);
            this->smList.push_back(sm);
            lineID = -2;
        }
    }
    input.close();

    for(int i=0;i<this->bsList.size();i++){
        this->involvedPairs.push_back(new vector<int>());
    }
    for(int i=0;i<this->bsPairs.size();i++){
        BackboneSitesPair* pair = bsPairs[i];
        int idA = pair->siteA->resseq;
        int idB = pair->siteB->resseq;
        this->involvedPairs[idA]->push_back(i);
        this->involvedPairs[idB]->push_back(i);
    }
}

void AADesignTemplate::loadS1S2(S1EnergyTable* s1Etable, S2MatrixFinder* s2Etable){
    int n = bsList.size();
    for(int i=0;i<n;i++){
        AAProbabilityArray* pa = new AAProbabilityArray();
        s1Etable->getS1(*bsList[i], pa);
        this->paList.push_back(pa);
    }

    int m = bsPairs.size();
    for(int i=0;i<m;i++){
        int seqSep = bsPairs[i]->seqSep;
        AAScoreMatrix* sm = new AAScoreMatrix();
        s2Etable->getSM(bsPairs[i], sm);
    //    float wt = 0.5;
    //    if(seqSep < 5) wt = 0.4;
        if(seqSep < 4)
            sm->multiply(0.5);
        else
            sm->multiply(0.8);
        this->smList.push_back(sm);
    }
    for(int i=0;i<n;i++){
        this->involvedPairs.push_back(new vector<int>());
    }
    for(int i=0;i<m;i++){
        BackboneSitesPair* pair = bsPairs[i];
        int idA = pair->siteA->resseq;
        int idB = pair->siteB->resseq;
        this->involvedPairs[idA]->push_back(i);
        this->involvedPairs[idB]->push_back(i);
    }
}

void AADesignTemplate::loadNearestPointS1S2(S1EnergyTable* s1Etable, S2MatrixFinder* s2Etable){
    int n = bsList.size();
    for(int i=0;i<n;i++){
        AAProbabilityArray* pa = new AAProbabilityArray();
        s1Etable->getS1NearestPoint(*bsList[i], pa);
        this->paList.push_back(pa);
    }

    int m = bsPairs.size();
    for(int i=0;i<m;i++){
        int seqSep = bsPairs[i]->seqSep;
        AAScoreMatrix* sm = new AAScoreMatrix();
        s2Etable->getSMNearestPair(bsPairs[i], sm);
        float wt = 0.5;
        if(seqSep < 5) wt = 0.4;
        sm->multiply(wt);
        this->smList.push_back(sm);
    }
    for(int i=0;i<n;i++){
        this->involvedPairs.push_back(new vector<int>());
    }
    for(int i=0;i<m;i++){
        BackboneSitesPair* pair = bsPairs[i];
        int idA = pair->siteA->resseq;
        int idB = pair->siteB->resseq;
        this->involvedPairs[idA]->push_back(i);
        this->involvedPairs[idB]->push_back(i);
    }
}

float AADesignTemplate::energy(AASequence* seq){
    if(seq->getLength() != this->paList.size()){
        cout << "sequence length not equal" << endl;
        exit(1);
    }
    int len = seq->getLength();
    double e = 0;
    for(int i=0;i<len;i++){
        double p = paList[i]->getProbability(seq->getChoice(i));
        e += -1.0*log(p);
    }

    for(int i=0;i<bsPairs.size();i++){
        BackboneSitesPair* bsPair = bsPairs[i];
        int idA = bsPair->siteA->resseq;
        int idB = bsPair->siteB->resseq;
        e += smList[i]->getValue(seq->getChoice(idA), seq->getChoice(idB));
    }
    return e;
}

float AADesignTemplate::mutEnergy(AASequence* oldSeq, int pos, int newChoice){
    int oldChoice = oldSeq->getChoice(pos);
    if(newChoice == oldChoice)
        return 0;
    float oldEnergy = -1.0*log(paList[pos]->getProbability(oldChoice));
    float newEnergy = -1.0*log(paList[pos]->getProbability(newChoice));

    vector<int>* pairIDs = this->involvedPairs[pos];
    for(int k=0;k<pairIDs->size();k++){
        int i = pairIDs->at(k);
        int idA = bsPairs.at(i)->siteA->resseq;
        int idB = bsPairs.at(i)->siteB->resseq;
        if(idA == pos){
            oldEnergy += smList[i]->getValue(oldChoice, oldSeq->getChoice(idB));
            newEnergy += smList[i]->getValue(newChoice, oldSeq->getChoice(idB));
        }
        else if(idB == pos){
            oldEnergy += smList[i]->getValue(oldSeq->getChoice(idA), oldChoice);
            newEnergy += smList[i]->getValue(oldSeq->getChoice(idA), newChoice);
        }
        else{
            cout << "IDA " << idA << " IDB " << idB << " pos " << pos << endl;
            exit(1);
        }
    }
    return newEnergy - oldEnergy;
}

float AADesignTemplate::mutEnergyS1(AASequence* oldSeq, int pos, int newChoice){
    int oldChoice = oldSeq->getChoice(pos);
    if(newChoice == oldChoice)
        return 0;
    float oldEnergy = -1.0*log(paList[pos]->getProbability(oldChoice));
    float newEnergy = -1.0*log(paList[pos]->getProbability(newChoice));
    return newEnergy - oldEnergy;
}

bool AADesignTemplate::accept(float deltaE, float T){
    double p = exp(-1.0*deltaE/T);
    double r = rand()/double(RAND_MAX);
    return r < p;
}

string AADesignTemplate::mcMinimize(){
    int len = paList.size();
    int stepNum = 10000;
    AASequence seq(paList.size());
    srand((unsigned)time(NULL));
    for(float T=10.0;T>0.0001;T=T*0.95){
        for(int i=0;i<stepNum;i++){
            int pos = rand()%len;
            int mutType = rand()%20;
            float me = mutEnergy(&seq, pos, mutType);
            if(accept(me, T))
                seq.applyMutation(pos,mutType);
        }
    }
    return seq.toString();
}

void AADesignTemplate::fixTemperatureDesign(float T, int seqNumber, vector<string>& results){
//    cout << "start aa design:" << endl;
    int len = paList.size();
    int stepNum = 10000;
    AASequence seq(paList.size());
    srand((unsigned)time(NULL));

    for(int k=0;k<seqNumber;k++){
        for(int i=0;i<stepNum;i++){
            int pos = rand()%len;
            int mutType = rand()%20;
            float me = mutEnergy(&seq, pos, mutType);
            if(accept(me, T))
                seq.applyMutation(pos,mutType);
        }
        results.push_back(seq.toString());
    }
}

vector<int> AADesignTemplate::getNativeRankList(){

    vector<int> ranks;
    int len = paList.size();
    ResName rn;
    char s[len+1];
    for(int i=0;i<len;i++){
        char c = rn.triToSin(this->bsList.at(i)->resname);
        s[i] = c;
    }
    s[len] = '\0';
    string natSeq = string(s);
    AASequence  seq(natSeq);

    for(int pos=0;pos<len;pos++){
        int rank = 1;
        for(int i=0;i<20;i++){
            float me = mutEnergy(&seq, pos, i);
            if(me < 0)
                rank++;
        }
        ranks.push_back(rank);
    }
    return ranks;
}

vector<double> AADesignTemplate::getNativePListS1(){
    vector<double> pList;
    int len = paList.size();
    ResName rn;
    for(int pos=0;pos<len;pos++){
        int aaType = rn.triToInt(this->bsList.at(pos)->resname);
        double p = this->paList.at(pos)->getProbability(aaType);
        pList.push_back(p);
    }
    return pList;
}

vector<int> AADesignTemplate::getNativeRankListS1(){
    vector<int> ranks;
    int len = paList.size();
    ResName rn;
    char s[len+1];
    for(int i=0;i<len;i++){
        char c = rn.triToSin(this->bsList.at(i)->resname);
        s[i] = c;
    }
    s[len] = '\0';
    string natSeq = string(s);
    AASequence  seq(natSeq);

    for(int pos=0;pos<len;pos++){
        int rank = 1;
        for(int i=0;i<20;i++){
            float me = mutEnergyS1(&seq, pos, i);
            if(me < 0)
                rank++;
        }
        ranks.push_back(rank);
    }
    return ranks;
}

void AADesignTemplate::getProfile(SeqProfile* prof){
    vector<string> seqs;
    fixTemperatureDesign(1.2, 1000, seqs);
    prof->initProfile(seqs);
}

void AADesignTemplate::getProfile(SeqProfile* prof, float T){
    vector<string> seqs;
    fixTemperatureDesign(T, 1000, seqs);
    prof->initProfile(seqs);
}

AADesignTemplate::~AADesignTemplate() {

    int n = this->paList.size();
    for(int i=0;i<n;i++){
        delete this->paList[i];
        delete this->involvedPairs[i];
    }

    int m = this->bsPairs.size();
    for(int i=0;i<m;i++){
        delete this->bsPairs[i];
        delete this->smList[i];
    }


    // TODO Auto-generated destructor stub
}

} /* namespace NSPdesignseq */

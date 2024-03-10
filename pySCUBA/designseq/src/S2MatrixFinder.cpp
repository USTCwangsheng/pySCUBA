/*
 * S2MatrixFinder.cpp
 *
 *  Created on: 2017Äê10ÔÂ26ÈÕ
 *      Author: notxp
 */

#include "designseq/S2MatrixFinder.h"

namespace NSPdesignseq {

S2MatrixFinder::S2MatrixFinder() {
    string repPath = NSPdataio::datapath() + "S2EnergyTable/repPoints2/";

    loadAllRepPoints(repPath);

    loadAllSaiPoints(repPath);

}

void S2MatrixFinder::loadRepPoints(string key, string fileName){

     ifstream file;
     file.open(fileName.c_str(), ios::in);
     if(! file.is_open()){
         cout << "fail to open file " << fileName << endl;
         exit(1);
     }
     string s;
     string substr;
     vector<ResPairOrientation> rpList;
     while(getline(file,s)){
         if(s.length() < 100) continue;
         substr = s.substr(8,119);
         ResPairOrientation rpo(substr);
         rpList.push_back(rpo);
     }
     repRPs[key] = rpList;
}

/*

void S2MatrixFinder::loadAllRepPoints(string& path){
    loadRepPoints("HH1", path + "HH1-100");
    loadRepPoints("HE1", path + "HE1-100");
    loadRepPoints("HC1", path + "HC1-200");
    loadRepPoints("EH1", path + "EH1-100");
    loadRepPoints("EE1", path + "EE1-100");
    loadRepPoints("EC1", path + "EC1-200");
    loadRepPoints("CH1", path + "CH1-200");
    loadRepPoints("CE1", path + "CE1-200");
    loadRepPoints("CC1", path + "CC1-200");

    loadRepPoints("HH2", path + "HH2-100");
    loadRepPoints("HE2", path + "HE2-100");
    loadRepPoints("HC2", path + "HC2-300");
    loadRepPoints("EH2", path + "EH2-100");
    loadRepPoints("EE2", path + "EE2-100");
    loadRepPoints("EC2", path + "EC2-400");
    loadRepPoints("CH2", path + "CH2-300");
    loadRepPoints("CE2", path + "CE2-400");
    loadRepPoints("CC2", path + "CC2-400");

    loadRepPoints("HH3", path + "HH3-100");
    loadRepPoints("HE3", path + "HE3-100");
    loadRepPoints("HC3", path + "HC3-400");
    loadRepPoints("EH3", path + "EH3-100");
    loadRepPoints("EE3", path + "EE3-200");
    loadRepPoints("EC3", path + "EC3-600");
    loadRepPoints("CH3", path + "CH3-400");
    loadRepPoints("CE3", path + "CE3-600");
    loadRepPoints("CC3", path + "CC3-600");

    loadRepPoints("HH4", path + "HH4-100");
    loadRepPoints("HE4", path + "HE4-200");
    loadRepPoints("HC4", path + "HC4-800");
    loadRepPoints("EH4", path + "EH4-200");
    loadRepPoints("EE4", path + "EE4-200");
    loadRepPoints("EC4", path + "EC4-800");
    loadRepPoints("CH4", path + "CH4-800");
    loadRepPoints("CE4", path + "CE4-800");
    loadRepPoints("CC4", path + "CC4-1000");

    loadRepPoints("HH5", path + "HH5-2000");
    loadRepPoints("HE5", path + "HE5-2000");
    loadRepPoints("HC5", path + "HC5-3000");
    loadRepPoints("EH5", path + "EH5-2000");
    loadRepPoints("EE5", path + "EE5-2000");
    loadRepPoints("EC5", path + "EC5-3000");
    loadRepPoints("CH5", path + "CH5-3000");
    loadRepPoints("CE5", path + "CE5-3000");
    loadRepPoints("CC5", path + "CC5-3000");
}

*/

void S2MatrixFinder::loadAllRepPoints(string& path) {
    loadRepPoints("HH1", path + "HH1-100");
    loadRepPoints("HE1", path + "HE1-100");
    loadRepPoints("HC1", path + "HC1-200");
    loadRepPoints("EH1", path + "EH1-100");
    loadRepPoints("EE1", path + "EE1-100");
    loadRepPoints("EC1", path + "EC1-200");
    loadRepPoints("CH1", path + "CH1-200");
    loadRepPoints("CE1", path + "CE1-200");
    loadRepPoints("CC1", path + "CC1-200");

    loadRepPoints("HH2", path + "HH2-100");
    loadRepPoints("HE2", path + "HE2-100");
    loadRepPoints("HC2", path + "HC2-400");
    loadRepPoints("EH2", path + "EH2-100");
    loadRepPoints("EE2", path + "EE2-400");
    loadRepPoints("EC2", path + "EC2-600");
    loadRepPoints("CH2", path + "CH2-400");
    loadRepPoints("CE2", path + "CE2-600");
    loadRepPoints("CC2", path + "CC2-800");

    loadRepPoints("HH3", path + "HH3-200");
    loadRepPoints("HE3", path + "HE3-100");
    loadRepPoints("HC3", path + "HC3-400");
    loadRepPoints("EH3", path + "EH3-100");
    loadRepPoints("EE3", path + "EE3-200");
    loadRepPoints("EC3", path + "EC3-600");
    loadRepPoints("CH3", path + "CH3-400");
    loadRepPoints("CE3", path + "CE3-600");
    loadRepPoints("CC3", path + "CC3-1000");

    loadRepPoints("HH4", path + "HH4-100");
    loadRepPoints("HE4", path + "HE4-200");
    loadRepPoints("HC4", path + "HC4-800");
    loadRepPoints("EH4", path + "EH4-200");
    loadRepPoints("EE4", path + "EE4-200");
    loadRepPoints("EC4", path + "EC4-800");
    loadRepPoints("CH4", path + "CH4-800");
    loadRepPoints("CE4", path + "CE4-800");
    loadRepPoints("CC4", path + "CC4-1500");

    loadRepPoints("HH5", path + "HH5-5000");
    loadRepPoints("HE5", path + "HE5-5000");
    loadRepPoints("HC5", path + "HC5-6000");
    loadRepPoints("EH5", path + "EH5-5000");
    loadRepPoints("EE5", path + "EE5-3000");
    loadRepPoints("EC5", path + "EC5-5000");
    loadRepPoints("CH5", path + "CH5-6000");
    loadRepPoints("CE5", path + "CE5-5000");
    loadRepPoints("CC5", path + "CC5-8000");
}

void S2MatrixFinder::loadSaiPoints(string key, string fileName){
    vector<SaiPair> saiList;
    ifstream file;
    file.open(fileName.c_str(), ios::in);
    if(! file.is_open()){
        cout << "fail to open file " << fileName << endl;
        exit(1);
    }
    string s;


    vector<string> spt;
    while(getline(file,s)){
        if(s.length() < 10) continue;
        splitString(s, " ", &spt);
        double x = atof(spt[0].c_str());
        double y = atof(spt[1].c_str());
        SaiPair sp(x,y);
        saiList.push_back(sp);
    }
    this->repSAIs[key] = saiList;
}

void S2MatrixFinder::loadAllSaiPoints(string& path){



    string hec = "HEC";
    char tmp[100];
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            for(int k=1;k<5;k++)
            {
                sprintf(tmp, "%c%c%d", hec.at(i), hec.at(j), k);
                string key = string(tmp);
                sprintf(tmp, "%spair%d-sai",path.c_str(),k);
                string fileName = string(tmp);
                loadSaiPoints(key, fileName);
            }


            sprintf(tmp, "%c%c%d", hec.at(i), hec.at(j), 5);
            string key = string(tmp);
            string fileName = path + key + "-sai";
            loadSaiPoints(key, fileName);
        }
    }
}

string S2MatrixFinder::findMatrixFile(BackboneSitesPair& bsPair){
    char ssA = bsPair.siteA->sscode;
    char ssB = bsPair.siteB->sscode;
    int saiIndex = findSaiIndex(bsPair) + 1;
    char tmp[100];
    sprintf(tmp, "S2-%c%c%d-%d-4000", ssA, ssB, bsPair.seqSep, saiIndex);
    string key = string(tmp);
    return NSPdataio::datapath() + "S2EnergyTable/loose2/" + key;
}

string S2MatrixFinder::findMatrixFile(BackboneSitesPair& bsPair, int saiIndex){
    char ssA = bsPair.siteA->sscode;
    char ssB = bsPair.siteB->sscode;
    char tmp[100];
    sprintf(tmp, "S2-%c%c%d-%d-4000", ssA, ssB, bsPair.seqSep, saiIndex+1);
    string key = string(tmp);
    return NSPdataio::datapath() + "S2EnergyTable/loose2/" + key;
}

int S2MatrixFinder::findSaiIndex(BackboneSitesPair& bsPair){
//    cout << "find sai index: " << endl;
    char ssA = bsPair.siteA->sscode;
    char ssB = bsPair.siteB->sscode;
    int seqSep = bsPair.seqSep;
    char tmp[20];
    sprintf(tmp, "%c%c%d", ssA, ssB, seqSep);
    string key = string(tmp);
    double saiA = bsPair.siteA->data_[3];
    double saiB = bsPair.siteB->data_[3];
    SaiPair sp(saiA, saiB);

//    cout << "sai key: " << key << endl;


    vector<SaiPair>& list = this->repSAIs[key];
//    cout << "list size: " << list.size() << endl;

    double minDist = 10;
    int minIndex = -1;
    for(int i=0;i<12;i++){
        double d = sp.distanceSquare(&(list.at(i)));
        if(d < minDist){
            minDist = d;
            minIndex = i;
        }
    }
    return minIndex;
}

map<int,double> S2MatrixFinder::findSaiIndexList(BackboneSitesPair& bsPair){
    char ssA = bsPair.siteA->sscode;
    char ssB = bsPair.siteB->sscode;
    int seqSep = bsPair.seqSep;
    char tmp[20];
    sprintf(tmp, "%c%c%d", ssA, ssB, seqSep);
    string key = string(tmp);
    double saiA = bsPair.siteA->data_[3];
    double saiB = bsPair.siteB->data_[3];
    SaiPair sp(saiA, saiB);

//    cout << "sai key: " << key << endl;

    int n = 3;
    int k;
    int indexList[n];
    double distList[n];
    for(int i=0;i<n;i++){
        indexList[i] = -1;
        distList[i] = 999999.9;
    }

    vector<SaiPair>& list = this->repSAIs[key];
//    cout << "list size: " << list.size() << endl;

    double minDist = 10;
    int minIndex = -1;
    for(int i=0;i<12;i++){
        double d = sp.distanceSquare(&(list.at(i)));
        if(d < minDist){
            minDist = d;
            minIndex = i;
        }
    }
    map<int,double> saiMap;
    for(int i=0;i<12;i++){
        double d = sp.distanceSquare(&(list.at(i)));
        if(d<0.0001) d = 0.0001;
        if(d < distList[n-1]) {
              distList[n-1] = d;
              indexList[n-1] = i;
        }
        else
              continue;

        for(int j=n-2;j>=0;j--) {
              if(distList[j+1] < distList[j]) {
              d = distList[j];
              distList[j] = distList[j+1];
              distList[j+1] = d;
              k = indexList[j];
              indexList[j] = indexList[j+1];
              indexList[j+1] = k;
        }
        else
              break;
        }
    }
    for(int i=0;i<n;i++){
        saiMap[indexList[i]] = distList[i];
    }
    return saiMap;
}

int S2MatrixFinder::findRPOIndex(BackboneSitesPair& bsPair){
    char ssA = bsPair.siteA->sscode;
    char ssB = bsPair.siteB->sscode;
    int seqSep = bsPair.seqSep;
    char tmp[20];
    sprintf(tmp, "%c%c%d", ssA, ssB, seqSep);
    string key = string(tmp);

    ResPairOrientation rpo = bsPair.ori;
    vector<ResPairOrientation>& list = this->repRPs[key];
    double minDist = 10000;
    int minIndex = -1;
    for(int i=0;i<list.size();i++){
        double d = rpo.rmsd(list.at(i));
        if(d < minDist){
            minDist = d;
            minIndex = i;
        }
    }
    return minIndex;
}

map<int,double> S2MatrixFinder::findRPOIndexList(BackboneSitesPair& bsPair){

    int sep = bsPair.seqSep;


    char ssA = bsPair.siteA->sscode;
    char ssB = bsPair.siteB->sscode;
    int seqSep = bsPair.seqSep;
    char tmp[20];
    sprintf(tmp, "%c%c%d", ssA, ssB, seqSep);
    string key = string(tmp);

    ResPairOrientation rpo = bsPair.ori;
    vector<ResPairOrientation>& list = this->repRPs[key];
    int listSize = list.size();

    int n = 3;
    if(sep > 3)
        n = 4;

    int k;
    int indexList[n];
    double distList[n];
    for(int i=0;i<n;i++){
        indexList[i] = -1;
        distList[i] = 999999.9;
    }

    double minDist = 10000;
    vector<pair<int,double>> pairList;
    for(int i=0;i<listSize;i++){
        double d = rpo.rmsd(list[i]);
        if(d < 0.00001)
            d = 0.00001;
        if(d < distList[n-1]) {
            distList[n-1] = d;
            indexList[n-1] = i;
        }
        else
            continue;

        for(int j=n-2;j>=0;j--) {
            if(distList[j+1] < distList[j]) {
                d = distList[j];
                distList[j] = distList[j+1];
                distList[j+1] = d;
                k = indexList[j];
                indexList[j] = indexList[j+1];
                indexList[j+1] = k;
            }
            else
                break;
        }

    }
    map<int,double> rpoMap;
    for(int i=0;i<n;i++){
        rpoMap[indexList[i]] = distList[i];
    }
    return rpoMap;
}

void S2MatrixFinder::getSMNearestPair(BackboneSitesPair* rp, AAScoreMatrix* outputSM){

    string fileName = findMatrixFile(*rp);
    int rpoIndex = findRPOIndex(*rp);
    if(rpoIndex < 0){
        cerr << "can't find rp " << endl;
        exit(1);
    }

    int start = rpoIndex*22+2;
    ifstream file;
    file.open(fileName.c_str(), ios::in);
    if(! file.is_open()){
           cout << "fail to open file " << fileName << endl;
           exit(1);
    }
    string s;

    for(int i=0;i<start;i++)
        getline(file, s);
    vector<string> spt;
    float matrix[20][20];
    for(int i=0;i<20;i++){
        getline(file,s);
        splitString(s," ",&spt);
        for(int k=0;k<20;k++){
            matrix[i][k] = atof(spt[k].c_str());
        }
    }
    file.close();
    outputSM->initValue(matrix);
}

void S2MatrixFinder::getSM(BackboneSitesPair* rp, AAScoreMatrix* outputSM){
    map<int,double> saiMap = findSaiIndexList(*rp);
    map<int,double> rpoMap = findRPOIndexList(*rp);
    map<int,double>::iterator it1,it2;

    double sigmaRPO = 1.5;
    double sigmaSAI = 1.0;
    if(rp->seqSep == 5)
        sigmaRPO = 2.0;

//    cout << "saiList: " << saiMap.size() << endl;
//    cout << "rpoList: " << rpoMap.size() << endl;

    int maxRPOIndex = 0;
    for(it1=rpoMap.begin();it1!=rpoMap.end();++it1){

        if(it1->first > maxRPOIndex)
            maxRPOIndex = it1->first;
    }


    int rpoSites[maxRPOIndex+1];
    double rpoWTs[maxRPOIndex+1];
    for(int i=0;i<maxRPOIndex+1;i++)
    {
        rpoSites[i]=0;
        rpoWTs[i] = 0;
    }
    for(it1=rpoMap.begin();it1!=rpoMap.end();++it1){
        rpoSites[it1->first] = 1;
        rpoWTs[it1->first] = pow(it1->second, -sigmaRPO);
    }

    double wtSum = 0;
    double wt;
    for(it1=saiMap.begin();it1!=saiMap.end();++it1){
        for(it2=rpoMap.begin();it2!=rpoMap.end();++it2){
            wt = pow(it2->second,-sigmaRPO) * pow(it1->second,-sigmaSAI);
            wtSum += wt;

        }
    }


    float matrix[20][20];
    for(int i=0;i<20;i++){
        for(int j=0;j<20;j++){
            matrix[i][j] = 0;
        }
    }

    ifstream file;
    double wtSai, wtRPO;
    vector<string> spt;
    for(it1=saiMap.begin();it1!=saiMap.end();++it1){
        int saiIndex = it1->first;
        string fileName = findMatrixFile(*rp,saiIndex);
        ifstream file;
        file.open(fileName.c_str(), ios::in);
        if(! file.is_open()){
               cout << "fail to open file " << fileName << endl;
               exit(1);
        }
        string s;

        wtSai = pow(it1->second, -sigmaSAI);
        int lineID = -1;
        while(getline(file,s)){
            lineID++;
            int pairID = lineID/22;
            int typeA = lineID - pairID*22 -2;
            if(typeA < 0)
                continue;
            if(pairID > maxRPOIndex)
                break;
            if(rpoSites[pairID] == 0)
                continue;
            wt = wtSai*rpoWTs[pairID]/wtSum;
            splitString(s," ",&spt);
            for(int k=0;k<20;k++){
                matrix[typeA][k] += atof(spt[k].c_str()) * wt;
            }
        }
        file.close();
    }
    outputSM->initValue(matrix);
}


S2MatrixFinder::~S2MatrixFinder() {
    // TODO Auto-generated destructor stub
}

} /* namespace NSPdesignseq */

/*
 * S1EnergyTable.cpp
 *
 *  Created on: 2017Äê10ÔÂ26ÈÕ
 *      Author: notxp
 */

#include "designseq/S1EnergyTable.h"


namespace NSPdesignseq {



S1EnergyTable::S1EnergyTable() {

    string path = NSPdataio::datapath();
    string s1EnergyFile = path+"S1EnergyTable6";
    ifstream file;
    file.open(s1EnergyFile.c_str(), ios::in);
    if(! file.is_open()){
        cout << "fail to open file " << s1EnergyFile << endl;
        exit(1);
    }
    string s;
    vector<string> spt;
    char ssType;
    double sai, phi, psi;
    int ssInt, saiInt, ppInt;
    getline(file,s);
    while(getline(file,s)){
        if(s.length() < 174) continue;
        vector<string> spt;
        splitString(s," ", &spt);
        ssType = s.at(0);
        sai = atof(spt[1].c_str());
        phi = atof(spt[2].c_str());
        psi = atof(spt[3].c_str());
        ssInt = ssToInt(ssType);
        saiInt = saiToInt(sai);
        Phipsi pp(phi, psi);
        ppInt = ppLib.phipsiToIndex(&pp);
        for(int k=0;k<20;k++){
            double p = atof(spt[k+5].c_str());
            if(p == 0) p = 0.00001;
            s1ETable[ssInt][saiInt][ppInt][k] = p;
        }

    }

    saiList.push_back(0.1);
    for(int i=0;i<17;i++){
        saiList.push_back(0.05*i + 0.2);
    }

    // TODO Auto-generated constructor stub
}

vector<pair<int,double>> S1EnergyTable::getSaiList(double sai){
    vector<pair<int,double>> list;
    double d;
    for(int i=0;i<18;i++){
        d = (sai-saiList[i])*(sai-saiList[i]);
        if(d == 0)
            d = 0.0001;
        if(d < 0.01)
        {
            pair<int,double> p(i,d);
            list.push_back(p);
        }
    }
    return list;
}


void S1EnergyTable::getS1NearestPoint(NSPproteinrep::BackBoneSite& bs, AAProbabilityArray* pa){
    int ss = ssToInt(bs.sscode);
    int sai = saiToInt(bs.data_[3]);
    Phipsi pp(bs.data_[0], bs.data_[1]);
    int ppIndex = ppLib.phipsiToIndex(&pp);
    pa->initProbability(s1ETable[ss][sai][ppIndex]);
}


void S1EnergyTable::getS1(NSPproteinrep::BackBoneSite& bs, AAProbabilityArray* pa){
    int ss = ssToInt(bs.sscode);
    double sai = bs.data_[3];
    Phipsi pp(bs.data_[0], bs.data_[1]);
    vector<pair<int,double>> ppList = ppLib.neighborPhipsiIndexList(&pp);
    vector<pair<int,double>> saiList = getSaiList(sai);

    double wtSai,wtPP, dSai, dPP, wt;
    double wtSum = 0;
    int saiIndex, ppIndex;
    for(int i=0;i<saiList.size();i++){
        dSai = saiList[i].second;
        wtSai = 1.0/dSai;

        for(int j=0;j<ppList.size();j++){
            dPP = ppList[j].second;
            wtPP = 1.0/dPP/dPP;
            wtSum += wtSai*wtPP;

        }
    }

    double p[20];
    for(int i=0;i<20;i++){
        p[i] = 0;
    }


    for(int i=0;i<saiList.size();i++){
        dSai = saiList[i].second;
        saiIndex = saiList[i].first;
        wtSai = 1.0/dSai;
        for(int j=0;j<ppList.size();j++){
            dPP = ppList[j].second;
            ppIndex = ppList[j].first;
            wtPP = 1.0/dPP/dPP;
            wt = wtSai*wtPP/wtSum;

            for(int k=0;k<20;k++){
                p[k] += this->s1ETable[ss][saiIndex][ppIndex][k] * wt;
            }
        }
    }
    pa->initProbability(p);

}

S1EnergyTable::~S1EnergyTable() {
    // TODO Auto-generated destructor stub
}

} /* namespace NSPdesignseq */

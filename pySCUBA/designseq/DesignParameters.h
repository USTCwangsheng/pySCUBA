/*
 * DesignParameters.h
 *
 *  Created on: 2017年12月2日
 *      Author: notxp
 */

#ifndef DESIGNSEQ_DESIGNPARAMETERS_H_
#define DESIGNSEQ_DESIGNPARAMETERS_H_

#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include "designseq/StringTool.h"
#include "dataio/datapaths.h"
#include "designseq/ResName.h"
#include <cassert>

namespace NSPdesignseq {
using namespace std;

class DesignParameters {
public:

    double wtS1;
    double wtS2Local;
    double wtS2Nonlocal;
    double wtS2EE;

    double wtCorePacking;
    double wtSurfPacking;
    double packSigma;
    double packSAI0;

    double nonHBondPolarWeight;
    double repWeight;

    double ref[20];

    double refH[20];
    double refE[20];
    double refC[20];

    std::map<std::string, double> refLocal;

    double wtRotEnergy;

    double vdwLamda1;
    double vdwLamda2;

    double hbLamda1;
    double hbLamda2;

    double wdNP_NP;
    double wdP_NP;
    double wdP_P;

    double wdHB;
    double hbAngleFactor;
    double aroWD;
    double nb1Rescale;

    char designOrPacking;

    DesignParameters();
    DesignParameters(int tag);
    DesignParameters(string& paraFile);
    void readRef(string& refFile){
        ifstream f;
        f.open(refFile.c_str(),ios::in);
        string s;
        vector<string> spt;
        ResName rn;
        while(getline(f,s)){
            splitString(s," ",&spt);
            int aa = rn.triToInt(spt[0]);
            float e = atof(spt[1].c_str());
            this->ref[aa] = e;
        }
        f.close();
    }

    void readHECRef(string& refFile){
        ifstream f;
        f.open(refFile.c_str(),ios::in);
        if (!f.good())
        {
        	cout << "Reading refFile " << refFile.c_str() << " failed." << endl;
        	exit(1);
        }
        string s;
        vector<string> spt;
        ResName rn;
        while(getline(f,s)){
            splitString(s," ",&spt);
            int aa = rn.triToInt(spt[0]);
            float eH = atof(spt[1].c_str());
            float eE = atof(spt[2].c_str());
            float eC = atof(spt[3].c_str());

            this->refH[aa] = eH;
            this->refE[aa] = eE;
            this->refC[aa] = eC;
        }
        f.close();
    }

    void readRefELocal(const std::string & refELocalFile) {
        assert(!refELocalFile.empty());
        refLocal.clear();
        ifstream f;
        f.open(refELocalFile.c_str(),ios::in);
        string s;
        vector<string> spt;
        while(getline(f,s)){
            splitString(s," ",&spt);
            float e[9];
            for (int i = 0; i < 9; i++) {
                e[i] = atof(spt[i+1].c_str());
            }
            refLocal[encodeLocal('H', 0.1, spt[0])] = e[0];
            refLocal[encodeLocal('H', 0.5, spt[0])] = e[1];
            refLocal[encodeLocal('H', 0.9, spt[0])] = e[2];
            refLocal[encodeLocal('E', 0.1, spt[0])] = e[3];
            refLocal[encodeLocal('E', 0.5, spt[0])] = e[4];
            refLocal[encodeLocal('E', 0.9, spt[0])] = e[5];
            refLocal[encodeLocal('C', 0.1, spt[0])] = e[6];
            refLocal[encodeLocal('C', 0.5, spt[0])] = e[7];
            refLocal[encodeLocal('C', 0.9, spt[0])] = e[8];
        }
        assert(!refLocal.empty());
        f.close();
    }

    static std::string encodeLocal(char ssType, double sai, const std::string & resName) {
        /*
         * definition
         * H: helix
         * E: beta strand
         * C: coil
         * B: buried
         * I: intermediate
         * S: surface
         * ----------
         *      B       I       S
         * H    <=0.33  <=0.66  <=1.0
         * E    <=0.33  <=0.50  <=1.0
         * C    <=0.33  <=0.66  <=1.0
         */
        char saiType = 'B';
        switch (ssType) {
            case 'H':
                if (sai <= 0.33) saiType = 'B';
                else if (sai > 0.66) saiType = 'S';
                else saiType = 'I';
                break;
            case 'E':
                if (sai <= 0.33) saiType = 'B';
                else if (sai > 0.5) saiType = 'S';
                else saiType = 'I';
                break;
            case 'C':
                if (sai <= 0.33) saiType = 'B';
                else if (sai > 0.66) saiType = 'S';
                else saiType = 'I';
                break;
            default:
                std::cerr << "Undefined SS type [" << ssType << "] "
                        "while encoding index for reference energy "
                        "calculation" << std::endl;
                exit(1);
        }
        return resName + ssType + saiType;
    }

    double getRefELocal(char ssType, double sai, const std::string & resName) const {
        std::string key = encodeLocal(ssType, sai, resName);
        if (refLocal.find(key) == refLocal.end()) {
            std::cerr << "Can't get RefELocal for " << key << std::endl;
            exit(1);
        }
        return refLocal.at(key);
    }

    double getVdwWeight(double sai){
        return wtSurfPacking + ( wtCorePacking - wtSurfPacking)/(1+exp((sai-packSAI0)/packSigma));
    }
    virtual ~DesignParameters();
};

} /* namespace NSPdesignseq */

#endif /* DESIGNSEQ_DESIGNPARAMETERS_H_ */

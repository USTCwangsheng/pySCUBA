/*
 * AtomicEnergyCalcular.h
 *
 *  Created on: 2017年12月2日
 *      Author: notxp
 */

#ifndef DESIGNSEQ_ATOMICENERGYCALCULAR_H_
#define DESIGNSEQ_ATOMICENERGYCALCULAR_H_


#include <iostream>
#include <map>
#include "designseq/DesignParameters.h"
#include "designseq/ProteinRep.h"
#include "geometry/xyz.h"
#include "designseq/Angles.h"

namespace NSPdesignseq {

using namespace std;



class AtomicEnergyCalcular {
public:

    float lamdaRep;
    float lamdaAtr;
    float wdNP_NP;
    float wdP_NP;
    float wdP_P;
    float wdAro;


    float wdHB;

    float lamdaHBRep;
    float lamdaHBAtr;
    float hbAngleFactor;
    float nb1Rescale;

    /*
     * non-hbond
     * first dimension: d0, from 3.0 to 4.0
     * second dimension: d^2, from 0 to 64
     */
    float vdwEnergyTable1[250][3200];
    float vdwEnergyTable2[250][12800];

    /*
     * hbond
     * first dimension: d0, from 2.5 to 3.5
     * second dimension: d^2, from 0 to 30
     */

    //float hbEnergyTable[100][3000];


    map<string, int> donorTypeToIndex;
    map<string, int> acceptorTypeToIndex;

    float donerAngleAve[54];
    float donerAngleSdLeft[54];
    float donerAngleSdRight[54];

    float acceptorAngleAve[54];
    float acceptorAngleSdLeft[54];
    float acceptorAngleSdRight[54];

    float hbWD[54];
    float distAve[54];
    float distSd[54];

    float vdw_rescale{ 1.0 };

public:

    void RescaleVDW(float v) { vdw_rescale = v; }
    DesignParameters* para;
    AtomicEnergyCalcular();
    AtomicEnergyCalcular(DesignParameters* para, double vdwres = 1.0);
    void initVdwEnergyTable(double vdwres);
    void initHBParameter();
    int donorHbNameToIndex(string s) { return donorTypeToIndex[s];};
    int acceptorHbNameToIndex(string s) {return acceptorTypeToIndex[s];};
    float getEnergy(float d, float d0, float wd, float rescale){
        float e;
        if(d < d0)
        {
            /*
             * added by xuyang, 2020.9.5
             *
             */
            double newmin = vdw_rescale * d0;
            if (d > newmin) {
                return wd / rescale;
            }
            else {
                float u = (d - newmin) / newmin / lamdaRep;
                return u * u + wd / rescale;
            }

            //float u = (d-d0)/d0/lamdaRep;
            //return u*u+wd/rescale;
        }
        else
        {
            float u = (d-d0)/d0/lamdaAtr;
            if(u>4)
                e = 0.0;
            else
                e = wd*exp(-u*u);
            return e/rescale;
        }
        return e;
    }

    float getEnergy2(float dd, float d0, float wd, float rescale){
        float d0Square = d0*d0;
        int ddIndex = (int)(dd*200);
        if(ddIndex >= 12800)
        {
            cerr << "invalid dd: " << dd << endl;
            exit(1);
        }
        int d0Index = (int)(d0*100-150);
        if(d0Index >= 250)
        {
            cerr << "invalid d0: " << d0 << endl;
            exit(1);
        }
        if(dd < d0Square)
        {
            return vdwEnergyTable1[d0Index][ddIndex] + wd/rescale;
        }
        else
        {
            return wd/rescale * vdwEnergyTable2[d0Index][ddIndex];
        }
    }

    float getHBEnergy(float d, float d0, double wd){
            float e;
            if(d < d0)
            {
                float u = (d-d0)/d0/lamdaHBRep;
                return u*u+wd;
            }
            else
            {
                float u = (d-d0)/d0/lamdaHBAtr;
                if(u>4)
                    e = 0.0;
                else
                    e = wd*expf(-u*u);
                return e;
            }
    }

    float getHBWd() {
        return this->wdHB;
    }

    float hbEnergy(PolarAtom* a, PolarAtom* b);
    float getAtomEnergy(XYZ* coordA, XYZ* dirA, AtomProperty* apA,  XYZ* coordB, XYZ* dirB, AtomProperty* apB);
    float getAtomEnergyPacking(XYZ* coordA, XYZ* dirA, AtomProperty* apA,  XYZ* coordB, XYZ* dirB, AtomProperty* apB);
    float getHBAtomEnergyPacking(XYZ* coordA, XYZ* dirA, AtomProperty* apA,  XYZ* coordB, XYZ* dirB, AtomProperty* apB);
    float getAtomEnergyNoHb(XYZ* coordA, XYZ* dirA, AtomProperty* apA,  XYZ* coordB, XYZ* dirB, AtomProperty* apB);
    virtual ~AtomicEnergyCalcular();
};



} /* namespace NSPdesignseq */

#endif /* DESIGNSEQ_ATOMICENERGYCALCULAR_H_ */

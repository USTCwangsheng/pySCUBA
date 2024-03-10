/*
 * AtomicEnergyCalcular.cpp
 *
 *  Created on: 2017年12月2日
 *      Author: notxp
 */

#include "designseq/AtomicEnergyCalcular.h"

namespace NSPdesignseq {

AtomicEnergyCalcular::AtomicEnergyCalcular() {
    // TODO Auto-generated constructor stub
    this->lamdaRep = 0.1;
    this->lamdaAtr = 0.15;
    this->wdNP_NP = -1.3;
    this->wdP_NP = -0.7;
    this->wdP_P = 0.0;
    this->wdAro = -1.5;

    this->wdHB = -0.2;
    this->lamdaHBRep = 0.1;
    this->lamdaHBAtr = 0.15;
    this->nb1Rescale = 1.5;

    initVdwEnergyTable(1.0);
    initHBParameter();
}

AtomicEnergyCalcular::AtomicEnergyCalcular(DesignParameters* para, double vdwres){
    this->lamdaRep = para->vdwLamda1;
    this->lamdaAtr = para->vdwLamda2;
    this->wdNP_NP = para->wdNP_NP;
    this->wdP_NP = para->wdP_NP;
    this->wdP_P = para->wdP_P;
    this->wdHB = para->wdHB;
    this->wdAro = para->aroWD;
    this->nb1Rescale = para->nb1Rescale;

    this->lamdaHBRep = para->hbLamda1;
    this->lamdaHBAtr = para->hbLamda2;
    this->hbAngleFactor = para->hbAngleFactor;

    this->para = para;

    if(abs(lamdaRep) < 0.001)
        lamdaRep = 0.001;
    if(abs(lamdaAtr) < 0.001)
        lamdaAtr = 0.001;
    if(abs(lamdaHBRep) < 0.001)
        lamdaHBRep = 0.001;
    if(abs(lamdaHBAtr) < 0.001)
        lamdaHBAtr = 0.001;

    initVdwEnergyTable(vdwres);
    initHBParameter();
}

void AtomicEnergyCalcular::initVdwEnergyTable(double vdwres)
{
    RescaleVDW(vdwres);
    int i,j;
    float dd,d,d0,u;
    for(i=0;i<250;i++){
        d0 = 1.5+i*0.01 + 0.005;

        for(j=0;j<3200;j++){
            dd = j*0.005 + 0.0025;
            d = sqrt(dd);
            /*
             * added by xuyang, 2020.9.5
             *
             */
            double newmin = vdw_rescale * d0;
            if (d > newmin) {
                u = 0.0;
            }
            else {
                u = (d - newmin) / newmin / lamdaRep;
            }
            //u = (d-d0)/d0/lamdaRep;
            vdwEnergyTable1[i][j] = u*u;
        }
    }

    for(i=0;i<250;i++){
        d0 = 1.5+i*0.01 + 0.005;
        for(j=0;j<12800;j++){
            dd = j*0.005 + 0.0025;
            d = sqrt(dd);
            u = (d-d0)/d0/lamdaAtr;
            vdwEnergyTable2[i][j] = exp(-u*u);
        }
    }

}

void AtomicEnergyCalcular::initHBParameter()
{
    ResName rn;
    for(int i=0;i<20;i++)
    {
        string triName = rn.intToTri(i);
        this->donorTypeToIndex[triName+"-N"] = 0;
        this->acceptorTypeToIndex[triName+"-O"] = 0;
    }

    this->donorTypeToIndex["BB-N"] = 0;
    this->donorTypeToIndex["HIS-ND1"] = 1;
    this->donorTypeToIndex["HIS-NE2"] = 1;
    this->donorTypeToIndex["ARG-NH1"] = 2;
    this->donorTypeToIndex["ARG-NH2"] = 2;
    this->donorTypeToIndex["ARG-NE"] = 3;
    this->donorTypeToIndex["LYS-NZ"] = 4;
    this->donorTypeToIndex["TRP-NE1"] = 5;
    this->donorTypeToIndex["THR-OG1"] = 6;
    this->donorTypeToIndex["SER-OG"] = 6;
    this->donorTypeToIndex["TYR-OH"] = 7;
    this->donorTypeToIndex["GLN-NE2"] = 8;
    this->donorTypeToIndex["ASN-ND2"] = 8;

    this->acceptorTypeToIndex["BB-O"] = 0;
    this->acceptorTypeToIndex["GLU-OE1"] = 1;
    this->acceptorTypeToIndex["GLU-OE2"] = 1;
    this->acceptorTypeToIndex["ASP-OD1"] = 1;
    this->acceptorTypeToIndex["ASP-OD2"] = 1;
    this->acceptorTypeToIndex["GLN-OE1"] = 2;
    this->acceptorTypeToIndex["ASN-OD1"] = 2;
    this->acceptorTypeToIndex["SER-OG"] = 3;
    this->acceptorTypeToIndex["THR-OG1"] = 3;
    this->acceptorTypeToIndex["TYR-OH"] = 4;
    this->acceptorTypeToIndex["HIS-ND1"] = 5;
    this->acceptorTypeToIndex["HIS-NE2"] = 5;

    hbWD[0] =  -1.000;
    hbWD[1] = -2.800;
    hbWD[2] = -2.500;
    hbWD[3] = -1.898;
    hbWD[4] = -1.623;
    hbWD[5] = -1.410;
    hbWD[6] = -0.710;
    hbWD[7] = -1.251;
    hbWD[8] = -0.989;
    hbWD[9] = -1.781;
    hbWD[10] = -1.153;
    hbWD[11] = -1.516;
    hbWD[12] = -1.645;
    hbWD[13] = -1.268;
    hbWD[14] = -0.815;
    hbWD[15] = -0.745;
    hbWD[16] = -0.488;
    hbWD[17] = -0.250;
    hbWD[18] = -1.151;
    hbWD[19] = -1.659;
    hbWD[20] = -0.970;
    hbWD[21] = -0.762;
    hbWD[22] = -0.444;
    hbWD[23] = 0.0;
    hbWD[24] = -1.404;
    hbWD[25] = -1.199;
    hbWD[26] = -1.016;
    hbWD[27] = -1.049;
    hbWD[28] = -0.793;
    hbWD[29] = -0.783;
    hbWD[30] = -1.896;
    hbWD[31] = -1.843;
    hbWD[32] = -1.403;
    hbWD[33] = -1.998;
    hbWD[34] = -1.602;
    hbWD[35] = -1.068;
    hbWD[36] = -0.500;
    hbWD[37] = -1.451;
    hbWD[38] = -0.941;
    hbWD[39] = -0.877;
    hbWD[40] = -2.094;
    hbWD[41] = -1.206;
    hbWD[42] = -1.795;
    hbWD[43] = -2.014;
    hbWD[44] = -1.183;
    hbWD[45] = -1.762;
    hbWD[46] = -1.051;
    hbWD[47] = -1.600;
    hbWD[48] = -1.086;
    hbWD[49] = -0.825;
    hbWD[50] = -1.014;
    hbWD[51] = -1.058;
    hbWD[52] = -0.887;
    hbWD[53] = -0.199;





    distAve[ 0] = 2.9300;
    distSd[ 0]  = 0.1286;
    donerAngleAve[ 0] = 165.00;
    donerAngleSdLeft[ 0] = 30.66;
    donerAngleSdRight[ 0] = 11.85;
    acceptorAngleAve[ 0] = 156.00;
    acceptorAngleSdLeft[ 0] = 38.63;
    acceptorAngleSdRight[ 0] = 11.53;
    distAve[ 1] = 2.8900;
    distSd[ 1]  = 0.1331;
    donerAngleAve[ 1] = 164.00;
    donerAngleSdLeft[ 1] = 18.88;
    donerAngleSdRight[ 1] = 13.83;
    acceptorAngleAve[ 1] = 130.00;
    acceptorAngleSdLeft[ 1] = 16.95;
    acceptorAngleSdRight[ 1] = 23.53;
    distAve[ 2] = 2.9000;
    distSd[ 2]  = 0.1262;
    donerAngleAve[ 2] = 166.00;
    donerAngleSdLeft[ 2] = 21.78;
    donerAngleSdRight[ 2] = 12.40;
    acceptorAngleAve[ 2] = 134.00;
    acceptorAngleSdLeft[ 2] = 13.72;
    acceptorAngleSdRight[ 2] = 25.12;
    distAve[ 3] = 3.0500;
    distSd[ 3]  = 0.1201;
    donerAngleAve[ 3] = 163.00;
    donerAngleSdLeft[ 3] = 21.57;
    donerAngleSdRight[ 3] = 10.65;
    acceptorAngleAve[ 3] = 140.00;
    acceptorAngleSdLeft[ 3] = 19.34;
    acceptorAngleSdRight[ 3] = 11.29;
    distAve[ 4] = 2.9800;
    distSd[ 4]  = 0.1173;
    donerAngleAve[ 4] = 168.00;
    donerAngleSdLeft[ 4] = 39.78;
    donerAngleSdRight[ 4] = 10.65;
    acceptorAngleAve[ 4] = 122.00;
    acceptorAngleSdLeft[ 4] = 15.09;
    acceptorAngleSdRight[ 4] = 30.69;
    distAve[ 5] = 3.0000;
    distSd[ 5]  = 0.1154;
    donerAngleAve[ 5] = 164.00;
    donerAngleSdLeft[ 5] = 19.02;
    donerAngleSdRight[ 5] = 11.77;
    acceptorAngleAve[ 5] = 166.00;
    acceptorAngleSdLeft[ 5] = 22.06;
    acceptorAngleSdRight[ 5] = 11.13;
    distAve[ 6] = 2.8100;
    distSd[ 6]  = 0.1527;
    donerAngleAve[ 6] = 166.00;
    donerAngleSdLeft[ 6] = 28.54;
    donerAngleSdRight[ 6] = 11.61;
    acceptorAngleAve[ 6] = 137.00;
    acceptorAngleSdLeft[ 6] = 16.26;
    acceptorAngleSdRight[ 6] = 24.80;
    distAve[ 7] = 2.7400;
    distSd[ 7]  = 0.1756;
    donerAngleAve[ 7] = 169.00;
    donerAngleSdLeft[ 7] = 22.73;
    donerAngleSdRight[ 7] = 9.54;
    acceptorAngleAve[ 7] = 120.00;
    acceptorAngleSdLeft[ 7] = 17.16;
    acceptorAngleSdRight[ 7] = 27.67;
    distAve[ 8] = 2.8000;
    distSd[ 8]  = 0.1605;
    donerAngleAve[ 8] = 167.00;
    donerAngleSdLeft[ 8] = 27.61;
    donerAngleSdRight[ 8] = 11.13;
    acceptorAngleAve[ 8] = 136.00;
    acceptorAngleSdLeft[ 8] = 19.50;
    acceptorAngleSdRight[ 8] = 21.94;
    distAve[ 9] = 2.7500;
    distSd[ 9]  = 0.1549;
    donerAngleAve[ 9] = 166.00;
    donerAngleSdLeft[ 9] = 23.33;
    donerAngleSdRight[ 9] = 11.29;
    acceptorAngleAve[ 9] = 113.00;
    acceptorAngleSdLeft[ 9] = 10.84;
    acceptorAngleSdRight[ 9] = 21.47;
    distAve[10] = 2.7100;
    distSd[10]  = 0.1531;
    donerAngleAve[10] = 166.00;
    donerAngleSdLeft[10] = 25.14;
    donerAngleSdRight[10] = 11.77;
    acceptorAngleAve[10] = 116.00;
    acceptorAngleSdLeft[10] = 9.09;
    acceptorAngleSdRight[10] = 12.24;
    distAve[11] = 3.0900;
    distSd[11]  = 0.1739;
    donerAngleAve[11] = 136.00;
    donerAngleSdLeft[11] = 17.55;
    donerAngleSdRight[11] = 40.70;
    acceptorAngleAve[11] = 137.00;
    acceptorAngleSdLeft[11] = 17.25;
    acceptorAngleSdRight[11] = 39.75;
    distAve[12] = 2.9000;
    distSd[12]  = 0.1533;
    donerAngleAve[12] = 131.00;
    donerAngleSdLeft[12] = 26.71;
    donerAngleSdRight[12] = 19.72;
    acceptorAngleAve[12] = 140.00;
    acceptorAngleSdLeft[12] = 20.10;
    acceptorAngleSdRight[12] = 23.53;
    distAve[13] = 2.8800;
    distSd[13]  = 0.1710;
    donerAngleAve[13] = 117.00;
    donerAngleSdLeft[13] = 16.59;
    donerAngleSdRight[13] = 20.51;
    acceptorAngleAve[13] = 117.00;
    acceptorAngleSdLeft[13] = 13.91;
    acceptorAngleSdRight[13] = 28.30;
    distAve[14] = 2.8900;
    distSd[14]  = 0.1682;
    donerAngleAve[14] = 131.00;
    donerAngleSdLeft[14] = 27.38;
    donerAngleSdRight[14] = 14.95;
    acceptorAngleAve[14] = 137.00;
    acceptorAngleSdLeft[14] = 20.54;
    acceptorAngleSdRight[14] = 24.17;
    distAve[15] = 2.9400;
    distSd[15]  = 0.1604;
    donerAngleAve[15] = 132.00;
    donerAngleSdLeft[15] = 25.13;
    donerAngleSdRight[15] = 17.49;
    acceptorAngleAve[15] = 119.00;
    acceptorAngleSdLeft[15] = 15.09;
    acceptorAngleSdRight[15] = 29.73;
    distAve[16] = 2.9500;
    distSd[16]  = 0.1534;
    donerAngleAve[16] = 132.00;
    donerAngleSdLeft[16] = 24.29;
    donerAngleSdRight[16] = 15.42;
    acceptorAngleAve[16] = 120.00;
    acceptorAngleSdLeft[16] = 16.65;
    acceptorAngleSdRight[16] = 27.03;
    distAve[17] = 3.0400;
    distSd[17]  = 0.1475;
    donerAngleAve[17] = 134.00;
    donerAngleSdLeft[17] = 28.85;
    donerAngleSdRight[17] = 12.08;
    acceptorAngleAve[17] = 155.00;
    acceptorAngleSdLeft[17] = 40.67;
    acceptorAngleSdRight[17] = 20.99;
    distAve[18] = 2.8600;
    distSd[18]  = 0.1392;
    donerAngleAve[18] = 164.00;
    donerAngleSdLeft[18] = 26.16;
    donerAngleSdRight[18] = 10.97;
    acceptorAngleAve[18] = 143.00;
    acceptorAngleSdLeft[18] = 17.45;
    acceptorAngleSdRight[18] = 25.12;
    distAve[19] = 2.8500;
    distSd[19]  = 0.1451;
    donerAngleAve[19] = 168.00;
    donerAngleSdLeft[19] = 17.37;
    donerAngleSdRight[19] = 9.70;
    acceptorAngleAve[19] = 116.00;
    acceptorAngleSdLeft[19] = 11.82;
    acceptorAngleSdRight[19] = 22.58;
    distAve[20] = 2.8400;
    distSd[20]  = 0.1478;
    donerAngleAve[20] = 163.00;
    donerAngleSdLeft[20] = 24.85;
    donerAngleSdRight[20] = 11.93;
    acceptorAngleAve[20] = 142.00;
    acceptorAngleSdLeft[20] = 19.96;
    acceptorAngleSdRight[20] = 24.96;
    distAve[21] = 2.9000;
    distSd[21]  = 0.1441;
    donerAngleAve[21] = 166.00;
    donerAngleSdLeft[21] = 28.65;
    donerAngleSdRight[21] = 11.29;
    acceptorAngleAve[21] = 128.00;
    acceptorAngleSdLeft[21] = 15.86;
    acceptorAngleSdRight[21] = 24.96;
    distAve[22] = 3.0800;
    distSd[22]  = 0.1311;
    donerAngleAve[22] = 165.00;
    donerAngleSdLeft[22] = 39.58;
    donerAngleSdRight[22] = 12.24;
    acceptorAngleAve[22] = 109.00;
    acceptorAngleSdLeft[22] = 12.00;
    acceptorAngleSdRight[22] = 38.48;
    distAve[23] = 2.9800;
    distSd[23]  = 0.1408;
    donerAngleAve[23] = 163.00;
    donerAngleSdLeft[23] = 38.14;
    donerAngleSdRight[23] = 15.26;
    acceptorAngleAve[23] = 165.00;
    acceptorAngleSdLeft[23] = 38.68;
    acceptorAngleSdRight[23] = 12.40;
    distAve[24] = 2.8200;
    distSd[24]  = 0.2222;
    donerAngleAve[24] = 102.00;
    donerAngleSdLeft[24] = 15.26;
    donerAngleSdRight[24] = 32.28;
    acceptorAngleAve[24] = 144.00;
    acceptorAngleSdLeft[24] = 21.63;
    acceptorAngleSdRight[24] = 23.06;
    distAve[25] = 2.7800;
    distSd[25]  = 0.2538;
    donerAngleAve[25] = 105.00;
    donerAngleSdLeft[25] = 14.86;
    donerAngleSdRight[25] = 24.80;
    acceptorAngleAve[25] = 122.00;
    acceptorAngleSdLeft[25] = 22.17;
    acceptorAngleSdRight[25] = 26.39;
    distAve[26] = 2.8100;
    distSd[26]  = 0.2351;
    donerAngleAve[26] = 103.00;
    donerAngleSdLeft[26] = 15.38;
    donerAngleSdRight[26] = 24.17;
    acceptorAngleAve[26] = 131.00;
    acceptorAngleSdLeft[26] = 17.69;
    acceptorAngleSdRight[26] = 27.67;
    distAve[27] = 2.8600;
    distSd[27]  = 0.2380;
    donerAngleAve[27] = 108.00;
    donerAngleSdLeft[27] = 18.89;
    donerAngleSdRight[27] = 38.16;
    acceptorAngleAve[27] = 131.00;
    acceptorAngleSdLeft[27] = 18.88;
    acceptorAngleSdRight[27] = 21.15;
    distAve[28] = 2.9400;
    distSd[28]  = 0.2343;
    donerAngleAve[28] = 105.00;
    donerAngleSdLeft[28] = 19.84;
    donerAngleSdRight[28] = 26.24;
    acceptorAngleAve[28] = 123.00;
    acceptorAngleSdLeft[28] = 18.27;
    acceptorAngleSdRight[28] = 24.49;
    distAve[29] = 2.8900;
    distSd[29]  = 0.2410;
    donerAngleAve[29] = 106.00;
    donerAngleSdLeft[29] = 15.77;
    donerAngleSdRight[29] = 19.08;
    acceptorAngleAve[29] = 167.00;
    acceptorAngleSdLeft[29] = 37.94;
    acceptorAngleSdRight[29] = 10.97;
    distAve[30] = 2.8900;
    distSd[30]  = 0.1075;
    donerAngleAve[30] = 164.00;
    donerAngleSdLeft[30] = 16.80;
    donerAngleSdRight[30] = 13.04;
    acceptorAngleAve[30] = 135.00;
    acceptorAngleSdLeft[30] = 16.35;
    acceptorAngleSdRight[30] = 24.01;
    distAve[31] = 2.8600;
    distSd[31]  = 0.1152;
    donerAngleAve[31] = 167.00;
    donerAngleSdLeft[31] = 13.23;
    donerAngleSdRight[31] = 11.13;
    acceptorAngleAve[31] = 120.00;
    acceptorAngleSdLeft[31] = 13.61;
    acceptorAngleSdRight[31] = 29.42;
    distAve[32] = 2.8800;
    distSd[32]  = 0.1159;
    donerAngleAve[32] = 166.00;
    donerAngleSdLeft[32] = 20.95;
    donerAngleSdRight[32] = 11.61;
    acceptorAngleAve[32] = 136.00;
    acceptorAngleSdLeft[32] = 19.73;
    acceptorAngleSdRight[32] = 20.67;
    distAve[33] = 2.9100;
    distSd[33]  = 0.1088;
    donerAngleAve[33] = 168.00;
    donerAngleSdLeft[33] = 20.56;
    donerAngleSdRight[33] = 10.81;
    acceptorAngleAve[33] = 120.00;
    acceptorAngleSdLeft[33] = 14.63;
    acceptorAngleSdRight[33] = 25.28;
    distAve[34] = 2.9300;
    distSd[34]  = 0.0865;
    donerAngleAve[34] = 165.00;
    donerAngleSdLeft[34] = 16.68;
    donerAngleSdRight[34] = 12.08;
    acceptorAngleAve[34] = 118.00;
    acceptorAngleSdLeft[34] = 16.04;
    acceptorAngleSdRight[34] = 21.94;
    distAve[35] = 2.9900;
    distSd[35]  = 0.1052;
    donerAngleAve[35] = 166.00;
    donerAngleSdLeft[35] = 32.12;
    donerAngleSdRight[35] = 10.81;
    acceptorAngleAve[35] = 157.00;
    acceptorAngleSdLeft[35] = 29.34;
    acceptorAngleSdRight[35] = 20.99;
    distAve[36] = 2.7300;
    distSd[36]  = 0.1547;
    donerAngleAve[36] = 103.00;
    donerAngleSdLeft[36] = 8.87;
    donerAngleSdRight[36] = 22.42;
    acceptorAngleAve[36] = 141.00;
    acceptorAngleSdLeft[36] = 17.58;
    acceptorAngleSdRight[36] = 15.26;
    distAve[37] = 2.6600;
    distSd[37]  = 0.1665;
    donerAngleAve[37] = 109.00;
    donerAngleSdLeft[37] = 10.18;
    donerAngleSdRight[37] = 16.06;
    acceptorAngleAve[37] = 120.00;
    acceptorAngleSdLeft[37] = 16.20;
    acceptorAngleSdRight[37] = 23.69;
    distAve[38] = 2.6900;
    distSd[38]  = 0.1734;
    donerAngleAve[38] = 110.00;
    donerAngleSdLeft[38] = 11.39;
    donerAngleSdRight[38] = 20.51;
    acceptorAngleAve[38] = 129.00;
    acceptorAngleSdLeft[38] = 17.75;
    acceptorAngleSdRight[38] = 24.17;
    distAve[39] = 2.7400;
    distSd[39]  = 0.1665;
    donerAngleAve[39] = 116.00;
    donerAngleSdLeft[39] = 13.39;
    donerAngleSdRight[39] = 29.42;
    acceptorAngleAve[39] = 118.00;
    acceptorAngleSdLeft[39] = 15.58;
    acceptorAngleSdRight[39] = 20.03;
    distAve[40] = 2.7100;
    distSd[40]  = 0.1362;
    donerAngleAve[40] = 115.00;
    donerAngleSdLeft[40] = 12.83;
    donerAngleSdRight[40] = 22.42;
    acceptorAngleAve[40] = 115.00;
    acceptorAngleSdLeft[40] = 8.82;
    acceptorAngleSdRight[40] = 14.95;
    distAve[41] = 2.7600;
    distSd[41]  = 0.1569;
    donerAngleAve[41] = 113.00;
    donerAngleSdLeft[41] = 10.73;
    donerAngleSdRight[41] = 20.35;
    acceptorAngleAve[41] = 165.00;
    acceptorAngleSdLeft[41] = 24.54;
    acceptorAngleSdRight[41] = 12.40;
    distAve[42] = 2.6700;
    distSd[42]  = 0.1275;
    donerAngleAve[42] = 113.00;
    donerAngleSdLeft[42] = 6.39;
    donerAngleSdRight[42] = 11.45;
    acceptorAngleAve[42] = 134.00;
    acceptorAngleSdLeft[42] = 15.98;
    acceptorAngleSdRight[42] = 20.51;
    distAve[43] = 2.6100;
    distSd[43]  = 0.1463;
    donerAngleAve[43] = 115.00;
    donerAngleSdLeft[43] = 6.13;
    donerAngleSdRight[43] = 10.65;
    acceptorAngleAve[43] = 121.00;
    acceptorAngleSdLeft[43] = 13.18;
    acceptorAngleSdRight[43] = 21.47;
    distAve[44] = 2.6600;
    distSd[44]  = 0.1565;
    donerAngleAve[44] = 114.00;
    donerAngleSdLeft[44] = 7.45;
    donerAngleSdRight[44] = 12.08;
    acceptorAngleAve[44] = 126.00;
    acceptorAngleSdLeft[44] = 13.90;
    acceptorAngleSdRight[44] = 25.12;
    distAve[45] = 2.7100;
    distSd[45]  = 0.1423;
    donerAngleAve[45] = 115.00;
    donerAngleSdLeft[45] = 8.85;
    donerAngleSdRight[45] = 14.31;
    acceptorAngleAve[45] = 115.00;
    acceptorAngleSdLeft[45] = 13.32;
    acceptorAngleSdRight[45] = 22.58;
    distAve[46] = 2.7100;
    distSd[46]  = 0.1401;
    donerAngleAve[46] = 115.00;
    donerAngleSdLeft[46] = 7.83;
    donerAngleSdRight[46] = 15.74;
    acceptorAngleAve[46] = 116.00;
    acceptorAngleSdLeft[46] = 8.95;
    acceptorAngleSdRight[46] = 15.11;
    distAve[47] = 2.7100;
    distSd[47]  = 0.1527;
    donerAngleAve[47] = 116.00;
    donerAngleSdLeft[47] = 9.97;
    donerAngleSdRight[47] = 12.24;
    acceptorAngleAve[47] = 166.00;
    acceptorAngleSdLeft[47] = 24.70;
    acceptorAngleSdRight[47] = 11.29;
    distAve[48] = 2.9300;
    distSd[48]  = 0.1362;
    donerAngleAve[48] = 117.00;
    donerAngleSdLeft[48] = 13.98;
    donerAngleSdRight[48] = 19.88;
    acceptorAngleAve[48] = 135.00;
    acceptorAngleSdLeft[48] = 17.35;
    acceptorAngleSdRight[48] = 21.78;
    distAve[49] = 2.9200;
    distSd[49]  = 0.1574;
    donerAngleAve[49] = 119.00;
    donerAngleSdLeft[49] = 10.86;
    donerAngleSdRight[49] = 20.67;
    acceptorAngleAve[49] = 120.00;
    acceptorAngleSdLeft[49] = 16.32;
    acceptorAngleSdRight[49] = 29.73;
    distAve[50] = 2.9300;
    distSd[50]  = 0.1479;
    donerAngleAve[50] = 116.00;
    donerAngleSdLeft[50] = 11.19;
    donerAngleSdRight[50] = 22.42;
    acceptorAngleAve[50] = 121.00;
    acceptorAngleSdLeft[50] = 13.57;
    acceptorAngleSdRight[50] = 33.23;
    distAve[51] = 2.9500;
    distSd[51]  = 0.1663;
    donerAngleAve[51] = 117.00;
    donerAngleSdLeft[51] = 14.88;
    donerAngleSdRight[51] = 24.01;
    acceptorAngleAve[51] = 122.00;
    acceptorAngleSdLeft[51] = 16.69;
    acceptorAngleSdRight[51] = 29.89;
    distAve[52] = 2.9800;
    distSd[52]  = 0.1777;
    donerAngleAve[52] = 120.00;
    donerAngleSdLeft[52] = 14.36;
    donerAngleSdRight[52] = 22.90;
    acceptorAngleAve[52] = 115.00;
    acceptorAngleSdLeft[52] = 13.90;
    acceptorAngleSdRight[52] = 24.65;
    distAve[53] = 3.0200;
    distSd[53]  = 0.1634;
    donerAngleAve[53] = 118.00;
    donerAngleSdLeft[53] = 16.96;
    donerAngleSdRight[53] = 25.44;
    acceptorAngleAve[53] = 161.00;
    acceptorAngleSdLeft[53] = 42.12;
    acceptorAngleSdRight[53] = 16.22;


}


/*
float AtomicEnergyCalcular::getEnergy(Tern* coordA, Tern* dirA, AtomProperty* apA,  Tern* coordB, Tern* dirB, AtomProperty* apB)
{
    if(!coordA->neighborTo(*coordB,6.0))
        return 0.0;
    float wd = 0, rescale;
    bool hbonded = false;
    int neighbor = apA->connectNum + apB->connectNum -1;
    if(neighbor  == 1)
        rescale = 1.5;
    else
        rescale = neighbor;

    if(!apA->isPolar && !apB->isPolar)
        wd = wdNP_NP;
    else if((apA->isHDonor && apB->isHAcceptor) || (apA->isHAcceptor && apB->isHDonor))
    {
        wd = wdHB;
        hbonded = true;
        return 0.0;
    }
    else if(apA->isPolar && apB->isPolar)
        wd = wdP_P;
    else if(apA->isPolar || apB->isPolar)
        wd = wdP_NP;

    Tern ab = *coordB - *coordA;
    float d = ab.length();

    float r1, r2, ang1, ang2, cosAngle;
    r1 = apA->vdwRadius;
    r2 = apB->vdwRadius;

    if(apA->isAromatic)
    {

        if(abs(dirA->length() - 1) > 0.0001)
            cerr << "invalid dir: " << dirA->t[0] << " " << dirA->t[1] << " " << dirA->t[2] << endl;
        cosAngle = ab*(*dirA)/d;
        if(cosAngle > 1)
            cosAngle = 1;
        if(cosAngle < -1)
            cosAngle = -1;
        if(cosAngle > 1 || cosAngle < -1)
        {
            cerr << "invalid cosAngle: " << cosAngle << endl;
            cerr << "dir: " << dirA->x() << " " << dirA->y() << " " << dirA->z() << endl;
        }
        ang1 = acos(cosAngle)*180/PI;


        r1 = apA->getSp2Radius(ang1);
    //    cout << "Ar1: " << r1 << " Ang1: " << ang1 << endl;
    }

    if(apB->isAromatic)
    {
        if(abs(dirB->length() - 1) > 0.01)
            cerr << "invalid dir: " << dirB->t[0] << " " << dirB->t[1] << " " << dirB->t[2] << endl;
        cosAngle = ab*(*dirB)/d;
        if(cosAngle > 1)
            cosAngle = 1;
        if(cosAngle < -1)
            cosAngle = -1;
        if(cosAngle > 1 || cosAngle < -1)
        {
            cerr << "invalid cosAngle: " << cosAngle << endl;
            cerr << "ab: " << ab.x() << " " << ab.y() << " " << ab.z() << endl;
            cerr << "dir: " << dirB->x() << " " << dirB->y() << " " << dirB->z() << endl;
        }
        ang2 = acos(cosAngle)*180/PI;
        r2 = apB->getSp2Radius(ang2);

    }

    if(hbonded)
    {
        r1 = apA->polarRadius;
        r2 = apB->polarRadius;
        return getEnergy(d,r1+r2,wd,1.0);
    }
    else
        return getEnergy(d,r1+r2,wd,rescale);
}
*/

float AtomicEnergyCalcular::getAtomEnergy(XYZ* coordA, XYZ* dirA, AtomProperty* apA,  XYZ* coordB, XYZ* dirB, AtomProperty* apB)
{
    float dd = coordA->squaredDistance(*coordB);
    float d = sqrt(dd);


    if(apA->isHDonor && apB->isHAcceptor){
        int idDoner = donorTypeToIndex[apA->atomUniqueName];
        int idAcceptor = acceptorTypeToIndex[apB->atomUniqueName];
        float hbDistance = distAve[idDoner*6+idAcceptor];
        return getHBEnergy(d, hbDistance, this->wdHB);
    }
    else if(apA->isHAcceptor && apB->isHDonor){
        int idDoner = donorTypeToIndex[apB->atomUniqueName];
        int idAcceptor = acceptorTypeToIndex[apA->atomUniqueName];
        float hbDistance = distAve[idDoner*6+idAcceptor];
        return getHBEnergy(d, hbDistance, this->wdHB);
    }

    if(dd > 60)
        return 0.0;

    float wd = 0, rescale;
    int neighbor = apA->connectNum + apB->connectNum -1;
    if(neighbor == 1)
        rescale = 1.5;
    else if(neighbor == 2)
        rescale = 3.0;
    else if(neighbor == 3)
        rescale = 4.0;
    else if(neighbor == 4)
        rescale = 5.0;
    else
        rescale = 6.0;


    if(!apA->isPolar && !apB->isPolar)
    {
        wd = wdNP_NP;
        if(apA->isAromatic || apB->isAromatic)
            wd = wdAro;
    }
    else if(apA->isPolar && apB->isPolar)
        wd = wdP_P;
    else
        wd = wdP_NP;

    XYZ ab = *coordB - *coordA;

    float r1, r2, cosAngleSquare, tmp;
    r1 = apA->vdwRadius;
    r2 = apB->vdwRadius;


    if(apA->isAromatic)
    {
        tmp = ab * (*dirA);
        cosAngleSquare = tmp*tmp/dd;
        r1 = apA->getSp2Radius2(cosAngleSquare);
    }

    if(apB->isAromatic)
    {
        tmp = ab*(*dirB);
        cosAngleSquare = tmp*tmp/dd;
        r2 = apB->getSp2Radius2(cosAngleSquare);
    }

    return getEnergy2(dd,r1+r2,wd,rescale);
}

float AtomicEnergyCalcular::getAtomEnergyPacking(XYZ* coordA, XYZ* dirA, AtomProperty* apA,  XYZ* coordB, XYZ* dirB, AtomProperty* apB)
{
    float dd = coordA->squaredDistance(*coordB);
    float d = sqrt(dd);


    if(apA->isHDonor && apB->isHAcceptor){
        int idDoner = donorTypeToIndex[apA->atomUniqueName];
        int idAcceptor = acceptorTypeToIndex[apB->atomUniqueName];
        float hbDistance = distAve[idDoner*6+idAcceptor];
        return getHBEnergy(d, hbDistance, this->hbWD[idDoner*6+idAcceptor]);
    }
    else if(apA->isHAcceptor && apB->isHDonor){
        int idDoner = donorTypeToIndex[apB->atomUniqueName];
        int idAcceptor = acceptorTypeToIndex[apA->atomUniqueName];
        float hbDistance = distAve[idDoner*6+idAcceptor];
        return getHBEnergy(d, hbDistance, this->hbWD[idDoner*6+idAcceptor]);
    }

    if(dd > 60)
        return 0.0;

    float wd = 0, rescale;
    int neighbor = apA->connectNum + apB->connectNum -1;
    if(neighbor == 1)
        rescale = 1.5;
    else if(neighbor == 2)
        rescale = 3.0;
    else if(neighbor == 3)
        rescale = 4.0;
    else if(neighbor == 4)
        rescale = 5.0;
    else
        rescale = 6.0;


    if(!apA->isPolar && !apB->isPolar)
    {
        wd = wdNP_NP;
        if(apA->isAromatic || apB->isAromatic)
            wd = wdAro;
    }
    else if(apA->isPolar && apB->isPolar)
        wd = wdP_P;
    else
        wd = wdP_NP;

    XYZ ab = *coordB - *coordA;

    float r1, r2, cosAngleSquare, tmp;
    r1 = apA->vdwRadius;
    r2 = apB->vdwRadius;


    if(apA->isAromatic)
    {
        tmp = ab * (*dirA);
        cosAngleSquare = tmp*tmp/dd;
        r1 = apA->getSp2Radius2(cosAngleSquare);
    }

    if(apB->isAromatic)
    {
        tmp = ab*(*dirB);
        cosAngleSquare = tmp*tmp/dd;
        r2 = apB->getSp2Radius2(cosAngleSquare);
    }

    return getEnergy2(dd,r1+r2,wd,rescale);
}

float AtomicEnergyCalcular::getHBAtomEnergyPacking(XYZ* coordA, XYZ* dirA, AtomProperty* apA,  XYZ* coordB, XYZ* dirB, AtomProperty* apB)
{
    float dd = coordA->squaredDistance(*coordB);
    float d = sqrt(dd);


    if(apA->isHDonor && apB->isHAcceptor){
        int idDoner = donorTypeToIndex[apA->atomUniqueName];
        int idAcceptor = acceptorTypeToIndex[apB->atomUniqueName];
        float hbDistance = distAve[idDoner*6+idAcceptor];
        return getHBEnergy(d, hbDistance, this->hbWD[idDoner*6+idAcceptor]);
    }
    else if(apA->isHAcceptor && apB->isHDonor){
        int idDoner = donorTypeToIndex[apB->atomUniqueName];
        int idAcceptor = acceptorTypeToIndex[apA->atomUniqueName];
        float hbDistance = distAve[idDoner*6+idAcceptor];
        return getHBEnergy(d, hbDistance, this->hbWD[idDoner*6+idAcceptor]);
    }

    return 0.0;
}

float AtomicEnergyCalcular::getAtomEnergyNoHb(XYZ* coordA, XYZ* dirA, AtomProperty* apA,  XYZ* coordB, XYZ* dirB, AtomProperty* apB)
{
    float dd = coordA->squaredDistance(*coordB);
    float d = sqrt(dd);

    if(apA->isHDonor && apB->isHAcceptor){
        return 0;
    }
    else if(apA->isHAcceptor && apB->isHDonor){
        return 0;
    }

    if(dd > 49)
        return 0.0;

    float wd = 0, rescale;
    int neighbor = apA->connectNum + apB->connectNum -1;
    if(neighbor == 1)
        rescale = 1.0;
    else
        rescale = neighbor/1.5;


    if(!apA->isPolar && !apB->isPolar)
    {
        wd = wdNP_NP;
        if(apA->isAromatic || apB->isAromatic)
            wd = wdAro;
    }
    else if(apA->isPolar && apB->isPolar)
        wd = wdP_P;
    else
        wd = wdP_NP;

    XYZ ab = *coordB - *coordA;

    float r1, r2, cosAngleSquare, tmp;
    r1 = apA->vdwRadius;
    r2 = apB->vdwRadius;


    if(apA->isAromatic)
    {
        tmp = ab * (*dirA);
        cosAngleSquare = tmp*tmp/dd;
        r1 = apA->getSp2Radius2(cosAngleSquare);
    }

    if(apB->isAromatic)
    {
        tmp = ab*(*dirB);
        cosAngleSquare = tmp*tmp/dd;
        r2 = apB->getSp2Radius2(cosAngleSquare);
    }

    return getEnergy2(dd,r1+r2,wd,rescale);
}

float AtomicEnergyCalcular::hbEnergy(PolarAtom* a, PolarAtom* b)
{
    float d = a->getCore().distance(b->getCore());
    if(d > 3.5)
        return 0;
    int donorType, accType;
    float ang1, ang2;
    if(a->isDonerAtom() && b->isAcceptorAtom())
    {
        donorType = donorTypeToIndex[a->getName()];
        accType = acceptorTypeToIndex[b->getName()];
        ang1 = angleX(a->getSupport(), a->getCore(), b->getCore());
        ang2 = angleX(a->getCore(),b->getCore(),b->getSupport());
    }
    else if(a->isAcceptorAtom() && b->isDonerAtom())
    {
        donorType = donorTypeToIndex[b->getName()];
        accType = acceptorTypeToIndex[a->getName()];
        ang1 = angleX(b->getSupport(), b->getCore(), a->getCore());
        ang2 = angleX(b->getCore(), a->getCore(), a->getSupport());
    }
    else
    {
        return 0.0;
    }

    int pairType = donorType*6+accType;


    if(pairType > 54)
    {
        cout << "invalid pairType: " << pairType << endl;
        return 0.0;
    }
    float d0Av = distAve[pairType];
    float ang1Av = donerAngleAve[pairType];
    float ang1SdL = donerAngleSdLeft[pairType];
    float ang1SdR = donerAngleSdRight[pairType];
    float ang2Av = acceptorAngleAve[pairType];
    float ang2SdL = acceptorAngleSdLeft[pairType];
    float ang2SdR = acceptorAngleSdRight[pairType];

    float wd = this->hbWD[pairType];

    if(wd == 0)
    {
        d0Av += 0.5;
    }


    float e = 0;
    float eDist = getHBEnergy(d,d0Av,wd);

    float eAng1 =0;
    float eAng2 =0;

    if(ang1 < ang1Av)
    {
        float u = (ang1Av-ang1)/ang1SdL;
        eAng1 = hbAngleFactor*u*u;
    }
    else if(ang1Av < 130)
    {
        float u = (ang1 - ang1Av)/ang1SdR;
        eAng1 = hbAngleFactor*u*u;
    }

    if(ang2 < ang2Av)
    {
        float u = (ang2Av-ang2)/ang2SdL;
        eAng2 = hbAngleFactor*u*u;
    }
    else if(ang2Av < 130)
    {
        float u = (ang2 - ang2Av)/ang2SdR;
        eAng2 = hbAngleFactor*u*u;
    }
    e = eDist + eAng1 + eAng2;
//    printf("pairType: %2d Dist: %5.3f Ang1: %7.3f Ang2: %7.3f eDist: %7.3f eAng1: %7.3f eAng2: %7.3f\n",pairType, d, ang1, ang2, eDist, eAng1, eAng2 );
    if(d > (a->getVdwRadius()+b->getVdwRadius()) && e > 0)
        return 0;
    return e;

}

AtomicEnergyCalcular::~AtomicEnergyCalcular() {
    // TODO Auto-generated destructor stub
}

} /* namespace xenergy */


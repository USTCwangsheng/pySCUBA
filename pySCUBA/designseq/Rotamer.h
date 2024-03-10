/*
 * Rotamer.h
 *
 *  Created on: 2017Äê10ÔÂ23ÈÕ
 *      Author: notxp
 */

#ifndef DESIGNSEQ_ROTAMER_H_
#define DESIGNSEQ_ROTAMER_H_

#include <string>
#include <vector>
#include "geometry/xyz.h"
#include "geometry/localframe.h"
#include "designseq/AtomLib.h"


namespace NSPdesignseq {
using namespace std;

class Rotamer {

public:
    string rotName;
    string triName;
    vector<string> atomNameList;
    vector<NSPgeometry::XYZ> coordList;
    bool hasAromaticRing;
    NSPgeometry::XYZ normalVectorOfRing;


    Rotamer();
    Rotamer(const string& rotName, const string& triName);

    void addAtom(string atomName, NSPgeometry::XYZ localCoord);
    NSPgeometry::XYZ& getAtomCoord(const string& atomName);
    void updateLawOfRing();
    void buildSidechain(NSPgeometry::LocalFrame& cs, vector<NSPgeometry::XYZ>& xyzList);
    double rms(Rotamer* other){
        if(this->coordList.size() != other->coordList.size())
            return 0;
        if(this->coordList.size() == 0)
            return 0;
        double tot = 0;
        for(int i=0;i<this->coordList.size();i++){
            tot += this->coordList[i].squaredDistance(other->coordList[i]);
        }
        return sqrt(tot/this->coordList.size());
    }
    virtual ~Rotamer();
};

} /* namespace NSPdesignseq */

#endif /* DESIGNSEQ_ROTAMER_H_ */

/*
 * Rotamer.cpp
 *
 *  Created on: 2017Äê10ÔÂ23ÈÕ
 *      Author: notxp
 */

#include "designseq/Rotamer.h"

namespace NSPdesignseq {

Rotamer::Rotamer() {
    // TODO Auto-generated constructor stub
    this->triName = "XXX";
    this->rotName = "XXX-0";
    this->hasAromaticRing = false;
}

Rotamer::Rotamer(const string& rotName, const string& triName)
{
    this->rotName = rotName;
    this->triName = triName;
    if(triName == "PHE" || triName == "TRP" || triName == "TYR" || triName == "HIS")
        this->hasAromaticRing = true;
    else
        this->hasAromaticRing = false;
}

void Rotamer::addAtom(string atomName, NSPgeometry::XYZ localCoord)
{
    this->atomNameList.push_back(atomName);
    this->coordList.push_back(localCoord);
}

NSPgeometry::XYZ& Rotamer::getAtomCoord(const string& atomName)
{
    for(unsigned int i=0;i<this->atomNameList.size();i++)
    {
        if(this->atomNameList.at(i) == atomName)
            return this->coordList.at(i);
    }
    cerr << "can't find atom: " << atomName << endl;
    exit(1);
}

void Rotamer::updateLawOfRing(){
    NSPgeometry::XYZ a,b,c;
    if(triName == "PHE" || triName == "TYR" || triName == "TRP")
    {
        a = getAtomCoord("CG");
        b = getAtomCoord("CD1");
        c = getAtomCoord("CD2");
        this->normalVectorOfRing = ~((b-a)^(c-a));
    }
    else if(triName == "HIS")
    {
        a = getAtomCoord("CG");
        b = getAtomCoord("ND1");
        c = getAtomCoord("CD2");
        this->normalVectorOfRing = ~((b-a)^(c-a));
    }
}

void Rotamer::buildSidechain(NSPgeometry::LocalFrame& cs, vector<NSPgeometry::XYZ>& xyzList){
    for(NSPgeometry::XYZ t : coordList){
        xyzList.push_back(cs.local2globalcrd(t));
    }
}

Rotamer::~Rotamer() {
    // TODO Auto-generated destructor stub
}

} /* namespace NSPdesignseq */

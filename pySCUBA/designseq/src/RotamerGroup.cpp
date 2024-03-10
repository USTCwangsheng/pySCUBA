/*
 * RotamerGroup.cpp
 *
 *  Created on: 2017年10月23日
 *      Author: notxp
 */

#include "designseq/RotamerGroup.h"

namespace NSPdesignseq {

RotamerGroup::RotamerGroup() {
    // TODO Auto-generated constructor stub
    for(int i=0;i<20;i++){
        aaRots.push_back(new vector<Rotamer*>());
    }
    this->rotNum = 0;
}


void RotamerGroup::addRotamer(Rotamer* rot)
{
    ResName rn;
    int aaType = rn.triToInt(rot->triName);
    this->rotList.push_back(rot);
    this->aaRots.at(aaType)->push_back(rot);
    this->rotNum ++;

    //char aaid = rn.intToSin(aaType);
    //if(idinlist.find(aaid)==idinlist.end()) idinlist.insert({aaid,std::set<int>()});
    //idinlist.at(aaid).insert(rotList.size()-1);
}

void RotamerGroup::deleteRotamer(int i){
    Rotamer* rot = rotList.at(i);
    string rotName = rot->rotName;
    ResName rn;
    int aaType = rn.triToInt(rot->triName);
    vector<Rotamer*>* aaRotList = aaRots.at(aaType);
    for(int k=0;k<aaRotList->size();k++){
        Rotamer* r = aaRotList->at(k);
        if(r->rotName == rot->rotName){
            aaRotList->erase(aaRotList->begin()+k);
            break;
        }
    }
    this->rotList.erase(rotList.begin()+i);
    this->rotNum --;

    //char aaid = rn.intToSin(aaType);
    //idinlist.at(aaid).erase(i);
    //for(auto &p:idinlist) {
    //    std::set<int> si;
    //    for(int ix:p.second) if(ix>i) si.insert(ix-1);
    //    p.second=si;
    //}
}

void RotamerGroup::merge(RotamerGroup& other, RotamerGroup& result){
    set<string> rotSet;
    for(int i=0;i<rotList.size();i++){
        rotSet.insert(rotList[i]->rotName);
        result.addRotamer(rotList[i]);
    }

    for(int i=0;i<other.rotList.size();i++){
        if(rotSet.find(other.rotList[i]->rotName) == rotSet.end()){
            result.addRotamer(other.rotList[i]);
        }
    }

}

void RotamerGroup::clear(){
    this->rotList.clear();
    this->rotNum = 0;
}

RotamerGroup::~RotamerGroup(){
    for(int i=0;i<20;i++){
        delete aaRots.at(i);
    }
}


} /* namespace NSPdesignseq */

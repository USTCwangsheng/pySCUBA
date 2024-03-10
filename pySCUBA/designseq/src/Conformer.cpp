/*
 * Conformer.cpp
 *
 *  Created on: 2018Äê1ÔÂ10ÈÕ
 *      Author: notxp
 */

#include "designseq/Conformer.h"

namespace NSPdesignseq {


void Conformer::addPolarAtoms(){
    int n = this->bbApList.size();
    string atomName;
    XYZ support, support1, support2;
    XYZ core;
    support = *this->bbCoordList[2];
    core = *this->bbCoordList[3];
    PolarAtom* pa = new PolarAtom(core,support,this->bbApList[3]);
    this->bbPolarList.push_back(pa);

    if(this->triName == "ASP" || this->triName == "ASN"){
        support = *this->scCoordList[1];
        core = *this->scCoordList[2];
        pa = new PolarAtom(core, support, this->scApList[2]);
        this->scPolarList.push_back(pa);
        core = *this->scCoordList[3];
        pa = new PolarAtom(core, support, this->scApList[2]);
        this->scPolarList.push_back(pa);
    }
    else if(this->triName == "GLU" || this->triName == "GLN"){
        support = *this->scCoordList[2];
        core = *this->scCoordList[3];
        pa = new PolarAtom(core, support, this->scApList[3]);
        this->scPolarList.push_back(pa);
        core = *this->scCoordList[4];
        pa = new PolarAtom(core, support, this->scApList[4]);
        this->scPolarList.push_back(pa);
    }
    else if(this->triName == "HIS"){
        core = *this->scCoordList[2];
        support1 = *this->scCoordList[1];
        support2 = *this->scCoordList[3];
        support = support1+support2-core;
        pa = new PolarAtom(core,support,this->scApList[2]);
        this->scPolarList.push_back(pa);

        core = *this->scCoordList[4];
        support1 = *this->scCoordList[3];
        support2 = *this->scCoordList[5];
        support = support1+support2-core;
        pa = new PolarAtom(core,support,this->scApList[4]);
        this->scPolarList.push_back(pa);
    }
    else if(this->triName == "LYS"){
        core = *this->scCoordList[4];
        support = *this->scCoordList[3];
        pa = new PolarAtom(core,support,this->scApList[4]);
        this->scPolarList.push_back(pa);
    }
    else if(this->triName == "ARG"){
        core = *this->scCoordList[3];
        support1 = *this->scCoordList[2];
        support2 = *this->scCoordList[4];
        support = support1+support2-core;
        pa = new PolarAtom(core,support,this->scApList[3]);
        this->scPolarList.push_back(pa);

        core = *this->scCoordList[5];
        support = *this->scCoordList[4];
        pa = new PolarAtom(core,support,this->scApList[5]);
        this->scPolarList.push_back(pa);

        core = *this->scCoordList[6];
        support = *this->scCoordList[4];
        pa = new PolarAtom(core,support,this->scApList[6]);
        this->scPolarList.push_back(pa);
    }
    else if(this->triName == "SER" || this->triName == "THR"){
        core = *this->scCoordList[1];
        support = *this->scCoordList[0];
        pa = new PolarAtom(core,support,this->scApList[1]);
        this->scPolarList.push_back(pa);
    }
    else if(this->triName == "TRP"){
        core = *this->scCoordList[3];
        support1 = *this->scCoordList[2];
        support2 = *this->scCoordList[4];
        support = support1 + support2 - core;
        pa = new PolarAtom(core,support,this->scApList[3]);
        this->scPolarList.push_back(pa);
    }
    else if(this->triName == "TYR"){
        core = *this->scCoordList[5];
        support = *this->scCoordList[4];
        pa = new PolarAtom(core,support,this->scApList[5]);
        this->scPolarList.push_back(pa);
    }
}

void Conformer::addPolarAtoms(BackBoneSite* bsPre){
    if(this->triName == "PRO"){
        return;
    }
    XYZ core = *this->bbCoordList[0];
    XYZ support1 = *this->bbCoordList[1];
    XYZ support2 = bsPre->ccrd();
    if(core.distance(support2) > 2) return;

    XYZ support = support1 + support2 - core;
    PolarAtom* pa = new PolarAtom(core, support, this->bbApList[0]);
    this->bbPolarList.push_back(pa);
}

Conformer::~Conformer(){
    XYZ* t;
    for(int i=0;i<this->bbCoordList.size();i++){
        t = this->bbCoordList.at(i);
        delete t;
    }
    for(int i=0;i<this->scCoordList.size();i++){
        t = this->scCoordList.at(i);
        delete t;
    }

    PolarAtom* p;
    for(int i=0;i<this->bbPolarList.size();i++){
        p = this->bbPolarList[i];
        delete p;
    }
    for(int i=0;i<this->scPolarList.size();i++){
        p = this->scPolarList[i];
        delete p;
    }
}


ConformerGroup::~ConformerGroup(){
    for(int i=0;i<this->confList.size();i++){
        delete confList.at(i);
    }
}
} /* namespace NSPdesignseq */

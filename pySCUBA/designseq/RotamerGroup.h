/*
 * RotamerGroup.h
 *
 *  Created on: 2017年10月23日
 *      Author: notxp
 */

#ifndef DESIGNSEQ_ROTAMERGROUP_H_
#define DESIGNSEQ_ROTAMERGROUP_H_

#include <vector>
#include <list>
#include <set>
#include "designseq/Rotamer.h"

namespace NSPdesignseq {

using namespace std;

class RotamerGroup {
public:
    vector<Rotamer*> rotList;
    vector<vector<Rotamer*>*> aaRots;

    /*
     * added by xuyang, 2019.4.5
     */
    std::map<char,std::set<int>> idinlist;



    int rotNum;
    RotamerGroup();
    void clear();
    void addRotamer(Rotamer* rot);
    void deleteRotamer(int i);

    void merge(RotamerGroup& other, RotamerGroup& result);


    /*
     * added by xuyang, 2019.4.5
     */
    bool containaa(char aaid) {
        if(idinlist.find(aaid)==idinlist.end()) return false;
        return true;
    }
    bool containother(char aaid) {
        for(auto &p:idinlist) {
            if(p.first==aaid) continue;
            return true;
        }
        return false;
    }
    void updataidinlist() {
        idinlist.clear();
        ResName rn;
        for(int i=0;i<rotList.size();i++) {
            char aaid = rn.triToSin(rotList[i]->triName);
            if(idinlist.find(aaid)==idinlist.end()) idinlist.insert({aaid,std::set<int>()});
            idinlist.at(aaid).insert(i);
        }
    }


    virtual ~RotamerGroup();
};

} /* namespace NSPdesignseq */

#endif /* DESIGNSEQ_ROTAMERGROUP_H_ */

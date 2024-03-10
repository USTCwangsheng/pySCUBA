/*
 * blkselector.cpp
 *
 *  Created on: 2019年12月12日
 *      Author: hyiu
 */
#include "iblock/blkselector.h"
#include "iblock/intrctmol.h"
#include  "dataio/inputlines.h"
#include <map>
#include <iostream>
using namespace NSPintrct;
using namespace NSPdataio;

void BlkSelector::addselection(const std::string & selectionstr){
	std::vector<std::string> words = NSPdataio::parseline(selectionstr,
			std::vector<int>());
	int cid = -1;
	std::map<int,std::string> chainselections;
	try {
		for (auto & w : words) {
			if (w.substr(0, 5) == "chain") {
				cid = std::stoi(w.substr(5));
				if (chainselections.find(cid) == chainselections.end()) {
					chainselections[cid] = std::string();
				}
			} else {
				assert(cid >= 0);
				chainselections.at(cid) += w + ",";
			}
		}
		for (auto &cs : chainselections) {
			std::vector<int> sel=NSPdataio::stringtoselection(cs.second);
			for(int i:sel)selectedblks[cs.first].insert(i);
		}
	} catch (std::exception &e) {
		std::cout << "Error procession selection string: " << selectionstr << std::endl;
		exit(1);
	}
	return;
}
std::vector<int> BlkSelector::selectatoms(const IntrctMol &imol) const{
	 std::vector<int> res;
	 for (auto &c:selectedblks){
		 for(auto r:c.second){
			 auto &blck=imol.getblck({c.first,r});
			 for(int i=blck.aoffset; i<blck.aoffset+blck.natoms();++i) res.push_back(i);
		 }
	 }
	return res;
}
std::vector<int> BlkSelector::selectmcatoms(const IntrctMol &imol) const{
	 std::vector<int> res;
	 if(selectedblks.empty()){
		 for(auto &bss:imol.bsinchains){
			 for(auto &bs:bss)
				 for(auto &i:bs.atomids())res.push_back(i);
		 }
		 return res;
	 }
	 for (auto &c:selectedblks){
		 for(auto r:c.second){
			 for(auto i:imol.getblck({c.first,r}).getbackboneatoms()) res.push_back(i);
		 }
	 }
	return res;
}
std::vector<int> BlkSelector::selectscatoms(const IntrctMol &imol) const{
	 std::vector<int> res;
	 if(selectedblks.empty()){
		 for(int c=0;c<imol.nchains();++c){
			 for(int r=0;r<imol.nresidues(c);++r){
				 auto & blck=imol.getblck({c,r});
				for(int i=blck.aoffset+2;i <blck.aoffset+blck.natoms()-2;++i) res.push_back(i);
			 }
		 }
		 return res;
	 }
	 for (auto &c:selectedblks){
		 for(auto r:c.second){
			 auto & blck=imol.getblck({c.first,r});
			 for(int i=blck.aoffset+2;i <blck.aoffset+blck.natoms()-2;++i) res.push_back(i);
		 }
	 }
	return res;
}

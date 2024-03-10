/*
 * isolatedpocket.cpp
 *
 *  Created on: 2020年7月13日
 *      Author: hyiu
 */
#include "proteinrep/isolatedpocket.h"
using namespace NSPproteinrep;
std::ostream & IsolatedPocket::writepdb(std::ostream &os){
	int atomid0=1;
	for (auto &r:components_){
		auto records=r.second.make_pdbrecords(r.first.chainid, r.first.residueid,
				' ', atomid0);
		atomid0 +=records.size();
		for(auto &ar:records) os<<ar.toString()<<std::endl;
	}
	return os;
}
void IsolatedPocket::addcomponents(const PdbReader &pdbreader,
		const std::vector<ComponentID> &cmpntids){
	AAConformersInModel pdbmodel;
	pdbmodel.getconformers(pdbreader);
	std::vector<ComponentID> cpts=cmpntids;
	if(cpts.empty()){
		for(auto &c: pdbmodel.conformers){
			for(auto &aa:c){
				cpts.push_back(ComponentID(aa));
			}
		}
	}
	for(auto &id:cpts){
		int chainnumber=pdbmodel.mappdbkeyint->chainNumber(id.chainid);
		int posinumber=pdbmodel.mappdbkeyint->posiNumber(std::make_pair(id.residueid,char(' ')),id.chainid);
		if(components_.find(id) != components_.end()){
			std::cout <<"Warning: component " << id.chainid <<" " << id.residueid <<" "<<id.resname <<
					"already exists" <<std::endl;
		}
		components_[id]=pdbmodel.conformers[chainnumber][posinumber];
		assert(id.chainid==components_[id].chainid_or);
		assert(id.residueid==components_[id].residueid_or);
		assert(id.resname==components_[id].residuename);
	}
}
IsolatedPocket::ComponentID IsolatedPocket::addcomponent(const AAConformer &aaconf,
		char newchainid,
		const std::vector<AtomID> & anchoratoms1,
		const std::vector<std::string> & anchoratomnms,
		const std::vector<double> &intcrds){
		int nrid=newresidueid(newchainid);
		ComponentID ncpid(newchainid,nrid,aaconf.residuename);
		components_[ncpid]=aaconf;
		auto & nconf=components_.at(ncpid);
		if(!anchoratoms1.empty()){
			assert(anchoratoms1.size()==3);
			assert(anchoratomnms.size()==3);
			for(auto &nm:anchoratomnms){
				assert(aaconf.getglobalcrd().find(nm) != aaconf.getglobalcrd().end());
			}
			std::vector<NSPgeometry::XYZ> joints;
			for(auto &a:anchoratoms1) joints.push_back(
					components_.at(a.cmpntID).getglobalcrd().at(a.atmnm));
			for(auto &a:anchoratomnms) joints.push_back(nconf.getglobalcrd().at(a));
			NSPgeometry::RigidTransform rt(joints,intcrds);
			for(auto &x:nconf.getglobalcrd()) rt.apply(&x.second);
		}
		return ncpid;
}


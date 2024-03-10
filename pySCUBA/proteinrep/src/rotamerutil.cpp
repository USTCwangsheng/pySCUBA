/*
 * rotamerutil.cpp
 *
 *  Created on: 2020年7月16日
 *      Author: hyiu
 */
#include "proteinrep/rotamerutil.h"
#include "pdbstatistics/phipsidistr.h"
#include "dstl/randomengine.h"
using namespace NSPdesignseq;
using namespace NSPproteinrep;
PhipsiLib & NSPproteinrep::getpplibinstance(){
	static PhipsiLib pplib;
	return pplib;
}
RotamerDistribution & RotamerDistribution::getdistribution(const std::string &residuetype){
	static std::map<std::string,RotamerDistribution> distributions;
	if(distributions.find(residuetype) == distributions.end()){
		distributions[residuetype]=RotamerDistribution();
		distributions.at(residuetype).setdistribution(residuetype);
	}
	return distributions.at(residuetype);
}
static std::vector<std::pair<double,double>> & phipsi200values(){
	static std::vector<std::pair<double,double>> values;
	if(values.empty()){
		auto &phipsidistr=NSPpdbstatistics::PhiPsiDistr::coildistr();
		double phi,psi;
		auto & rng=NSPdstl::RandomEngine<>::getinstance().realrng(0,1);
		for(int i=0;i<200;i++){
			phipsidistr.randomphipsi(rng,&phi,&psi);
			values.push_back(std::make_pair(phi,psi));
		}
	}
	return values;
}

void RotamerDistribution::setdistribution(const std::string &residuetype){
	rotamers.clear();
	rotamerene.clear();
	rotamerprob.clear();
	auto & pplib=getpplibinstance();
	std::string res=residuetype;
	RotamerLib & libA1 = RotamerLib::rotamerlib("A1");
	RotamerLib & libB1 = RotamerLib::rotamerlib("B1");
	RotamerGroup *aagrpA1 = libA1.getAAGroup(res);
	RotamerGroup *aagrpB1 = libB1.getAAGroup(res);
	auto & phipsivalues=phipsi200values();
	for(auto &r:aagrpA1->rotList){
		double ene=0.0;
		int count=0;
		for(auto &v:phipsivalues){
			Phipsi pp(v.first,v.second);
			if(pp.regionAB() !='A') continue;
			count++;
			ene+=libA1.getRotamerEnergy(r->rotName, pplib.phipsiToIndex(&pp));
		}
		ene=ene/(double)count;
		rotamerene.push_back(ene);
		rotamerprob.push_back(exp(-ene));
	}
	for(auto &r:aagrpB1->rotList){
		double ene=0.0;
		int count=0;
		for(auto &v:phipsivalues){
			Phipsi pp(v.first,v.second);
			if(pp.regionAB() !='B') continue;
			count++;
			ene+=libB1.getRotamerEnergy(r->rotName, pplib.phipsiToIndex(&pp));
		}
		rotamers.push_back(r);
		ene=ene/(double)count;
		rotamerene.push_back(ene);
		rotamerprob.push_back(exp(-ene));
	}
	double w=0.0;
	for(auto p:rotamerprob) w+=p;
	for(auto &p:rotamerprob) p=p/w;
}
void RotamerDistribution::setdistribution(const std::string &residuetype,double phi,double psi){
	rotamers.clear();
	rotamerene.clear();
	rotamerprob.clear();
    std::string res=residuetype;
    auto &pplib=getpplibinstance();
    Phipsi pp(phi,psi);
	RotamerLib & lib=RotamerLib::rotamerlibpp(phi,psi);
	RotamerGroup *aagrp = lib.getAAGroup(res);
	std::vector<double> prob;
	for (auto r : aagrp->rotList) {
		std::string rname = r->rotName;
		double rene = lib.getRotamerEnergy(rname, pplib.phipsiToIndex(&pp));
		rotamers.push_back(r);
		rotamerene.push_back(rene);
		rotamerprob.push_back(exp(-rene));
	}
	double w=0.0;
	for(auto p:rotamerprob) w+=p;
	for(auto &p:rotamerprob) p=p/w;
}

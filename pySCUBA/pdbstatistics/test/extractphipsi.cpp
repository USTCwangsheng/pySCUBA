/*
 * extractphipsi.cpp
 *
 *  Created on: 2017年1月17日
 *      Author: hyliu
 */

#include <backbone/backbonesite.h>
#include <iostream>
#include <fstream>
using namespace NSPproteinrep;
int main(int argc, char **argv) {
	std::vector<BackBoneSite> sites;
	readbackbonesites(std::string(argv[1]),sites);
	std::vector<std::pair<double,double>> phipsi_vals;
	std::vector<std::pair<double,double>> gly_phipsi;
	std::vector<std::pair<double,double>> prepro_phipsi;
	std::vector<std::pair<double,double>> transpro_phipsi;
	std::vector<std::pair<double,double>> cispro_phipsi;
	std::vector<std::pair<double,double>> helix_phipsi;
	std::vector<std::pair<double,double>> strand_phipsi;
	bool iscis=false;
	for(auto iter=sites.begin();iter != sites.end(); ++iter) {
		bool isgly=false;
		bool isprepro=false;
		bool istranspro=false;
		bool iscispro=false;
		bool ishelix=false;
		bool isstrand=false;
		bool iscoil=true;
		if(iter->sscode=='H') {
			ishelix=true;
			iscoil=false;
		} else if (iter->sscode=='E'){
			isstrand=true;
			iscoil=false;
		}
		if(iter->resname=="GLY") isgly=true;
		if(iter->resname=="PRO") {
			if(iscis)
				iscispro=true;
			else
				istranspro=true;
		}
		iscis=false;
		if( !chainendsite(iter,sites.end())) {
			auto next=iter+1;
			if(next->resname=="PRO") isprepro=true;
			if(iter->omiga() >-90 && iter->omiga() <90) iscis=true;
		}
		double phi=iter->phi();
		if(phi == 360.0) continue;
		double psi=iter->psi();
		if(psi ==360.0) continue;
		if(isgly) {
			gly_phipsi.push_back(std::make_pair(phi,psi));
		} else if(isprepro) {
			prepro_phipsi.push_back(std::make_pair(phi,psi));
		} else if(istranspro){
			transpro_phipsi.push_back(std::make_pair(phi,psi));
		} else if (iscispro){
			cispro_phipsi.push_back(std::make_pair(phi,psi));
		}
		else{
			if(ishelix) helix_phipsi.push_back(std::make_pair(phi,psi));
			else if (isstrand) strand_phipsi.push_back(std::make_pair(phi,psi));
			else phipsi_vals.push_back(std::make_pair(phi,psi));
		}
	}


	std::cout <<helix_phipsi.size()<<":"<<strand_phipsi.size()<<":"<<phipsi_vals.size() <<":"<<gly_phipsi.size()<<":" <<prepro_phipsi.size()
			<<":" <<transpro_phipsi.size()<<":"<<cispro_phipsi.size()<<std::endl;
	std::ofstream ofs;
	ofs.open("helixtorsions.dat");
	for(auto &pp:helix_phipsi){
		ofs <<pp.first <<"\t"<<pp.second<<std::endl;
	}
	ofs.close();
	ofs.open("strandtorsions.dat");
	for(auto &pp:strand_phipsi){
		ofs <<pp.first <<"\t"<<pp.second<<std::endl;
	}
	ofs.close();
	ofs.open("coiltorsions.dat");
	for(auto &pp:phipsi_vals){
		ofs <<pp.first <<"\t"<<pp.second<<std::endl;
	}
	ofs.close();
	ofs.open("glytorsions.dat");
	for(auto &pp:gly_phipsi){
		ofs <<pp.first <<"\t"<<pp.second<<std::endl;
	}
	ofs.close();
	ofs.open("preprotorsions.dat");
	for(auto &pp:prepro_phipsi){
		ofs <<pp.first <<"\t"<<pp.second<<std::endl;
	}
	ofs.close();
	ofs.open("transprotorsions.dat");
	for(auto &pp:transpro_phipsi){
		ofs <<pp.first <<"\t"<<pp.second<<std::endl;
	}
	ofs.close();
	ofs.open("cisprotorsions.dat");
	for(auto &pp:cispro_phipsi){
		ofs <<pp.first <<"\t"<<pp.second<<std::endl;
	}
	ofs.close();
}



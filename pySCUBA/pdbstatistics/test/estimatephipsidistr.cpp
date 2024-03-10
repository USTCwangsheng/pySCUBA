/*
 * estimatephipsidistr.cpp
 *
 *  Created on: 2016年11月27日
 *      Author: hyliu
 */
#include <backbone/backbonesite.h>
#include <pdbstatistics/phipsidistr.h>
#include <pdbstatistics/proteinblock.h>
#include <iostream>
#include <fstream>
using namespace NSPpdbstatistics;
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
		double phi=iter->phi();
		if(phi == 360.0) continue;
		double psi=iter->psi();
		if(psi ==360.0) continue;
		char pbtype=ProteinBlock::pbtype(sites,iter-sites.begin());
		if(pbtype=='m' && iter->sscode=='H'){
			helix_phipsi.push_back(std::make_pair(phi,psi));
			continue;
		}
		if(pbtype=='d' && iter->sscode=='E'){
			strand_phipsi.push_back(std::make_pair(phi,psi));
			continue;
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
			phipsi_vals.push_back(std::make_pair(phi,psi));
		}
	}
	PhiPsiDistr dist1(PhiPsiDistr::coilstep);
	PhiPsiDistr distgly(PhiPsiDistr::glystep);
	PhiPsiDistr distprepro(PhiPsiDistr::preprostep);
	PhiPsiDistr disttranspro(PhiPsiDistr::transprostep);
	PhiPsiDistr distcispro(PhiPsiDistr::cisprostep);
	PhiPsiDistr disthelix(PhiPsiDistr::coilstep);
	PhiPsiDistr diststrand(PhiPsiDistr::coilstep);
	std::cout <<phipsi_vals.size() <<":"<<gly_phipsi.size()<<":" <<prepro_phipsi.size()
			<<":" <<transpro_phipsi.size()<<":"<<cispro_phipsi.size()<<std::endl;
	std::cout <<helix_phipsi.size() <<":"<<strand_phipsi.size()<<std::endl;

	std::ofstream ofs;
	dist1.estimatedistr(phipsi_vals);
	ofs.open("coilphipsi.dat");
	dist1.writedistr(ofs);
	ofs.close();
	ofs.open("glyphipsi.dat");
	distgly.estimatedistr(gly_phipsi);
	distgly.writedistr(ofs);
	ofs.close();
	ofs.open("preprophipsi.dat");
	distprepro.estimatedistr(prepro_phipsi);
	distprepro.writedistr(ofs);
	ofs.close();
	ofs.open("transprophipsi.dat");
	disttranspro.estimatedistr(transpro_phipsi);
	disttranspro.writedistr(ofs);
	ofs.close();
	ofs.open("cisprophipsi.dat");
	distcispro.estimatedistr(cispro_phipsi);
	distcispro.writedistr(ofs);
	ofs.close();
	ofs.open("helixphipsi.dat");
	disthelix.estimatedistr(helix_phipsi);
	disthelix.writedistr(ofs);
	ofs.close();
	ofs.open("strandphipsi.dat");
	diststrand.estimatedistr(strand_phipsi);
	diststrand.writedistr(ofs);
	ofs.close();
	/*
	const PhiPsiDistr & dist1=PhiPsiDistr::coildistr();
	const PhiPsiDistr & distgly=PhiPsiDistr::glydistr();
	const PhiPsiDistr & distpro=PhiPsiDistr::preprodistr();
	double eav=0.0;
	double eav2=0.0;
	for(auto ps:phipsi_vals) {
		double e=dist1.statisticalenergy(ps.first,ps.second);
		eav += e;
		eav2 += e*e;
	}
	eav=eav/(double) (phipsi_vals.size());
	std::cout << eav <<" "<< sqrt(eav2/(double) (phipsi_vals.size())-eav*eav) <<std::endl;
	eav=0.0;
	eav2=0.0;
	for(auto ps:gly_phipsi) {
		double e=distgly.statisticalenergy(ps.first,ps.second);
		eav += e;
		eav2 += e*e;
	}
	eav=eav/(double) (gly_phipsi.size());
	std::cout << eav << " "<<sqrt(eav2/(double) (gly_phipsi.size())-eav*eav)<<std::endl;
	eav=0.0;
	eav2=0.0;
	for(auto ps:prepro_phipsi) {
		double e=distpro.statisticalenergy(ps.first,ps.second);
		eav += e;
		eav2 += e*e;
	}
	eav=eav/(double) (prepro_phipsi.size());
	std::cout << eav <<" "<< sqrt(eav2/(double) (prepro_phipsi.size())-eav*eav)<<std::endl;
	*/
}



/*
 * testphipsidistr.cpp
 *
 *  Created on: 2017年3月29日
 *      Author: hyliu
 */


#include "pdbstatistics/phipsidistr.h"
#include "backbone/backbonesite.h"
#include "dstl/randomengine.h"

using namespace NSPpdbstatistics;
using namespace NSPproteinrep;
using namespace NSPdstl;

int main(int argc, char **argv) {
	const PhiPsiDistr *coildistr=&(PhiPsiDistr::coildistr());
	const PhiPsiDistr *glydistr=&(PhiPsiDistr::glydistr());
	const PhiPsiDistr *preprodistr=&(PhiPsiDistr::preprodistr());
	const PhiPsiDistr *cisprodistr=&(PhiPsiDistr::cisprodistr());
	const PhiPsiDistr *transprodistr=&(PhiPsiDistr::transprodistr());
	const PhiPsiDistr *mixdistr=&(PhiPsiDistr::mixcoildistr());
	double step=0.1;
	for(int i=0;i<360/step;++i){
/*		for (int j=0;j<360/step;++j) {
			double phi=-180.0+step*(double) i;
			double psi=-180.0+step*(double) j;
			double ene=mixdistr->statisticalenergy(phi,psi);
			std::cout <<phi <<" "<< psi<<"  "<<ene <<std::endl;
		}*/
		double phi=step*(double) i;
		double dedphi;
		std::cout<< phi<<" "<< mixdistr->intplene_phi(phi,&dedphi) <<" "<<
				mixdistr->intplene_psi(phi,&dedphi)<<std::endl;
	}
//	coildistr->writedistr("phipsicoil.dat",-1.0);
/*	std::vector<BackBoneSite> sites;
	readbackbonesites(std::string(argv[1]),sites);
	RandomEngine<>::real_generator_type & rng= RandomEngine<>::getinstance().realrng();
	for(auto iter=sites.begin(); iter != sites.end();++iter) {
		double phi=iter->phi();
		double psi=iter->psi();

		if(phi ==360.0 || psi==360.0) continue;
		BackBoneSite * prev=nullptr;
		BackBoneSite *next=nullptr;
		if(! chainstartsite(iter) ) prev= &(*(iter-1));
		if(! chainendsite(iter,sites.end())) next=&(*(iter+1));
		const PhiPsiDistr *distr;
		std::string type;
		if(iter->resname == "GLY"){
			distr=glydistr;
			type="gly   ";
		} else if (iter->resname == "PRO") {
			if(prev &&(prev->omiga() > -90.0 and prev->omiga()<90.0)){
				type="cispro";
				distr=cisprodistr;
			}
			else {
				type="trspro";
				distr=transprodistr;
			}
		} else {
			if( next && next->resname == "PRO") {
				type="prepro";
				distr=preprodistr;
			}
			else {
				type="coil  ";
				distr=coildistr;
			}
		}
		double dedp,dedpsi;
		std::cout <<type<<"\t"<<phi<<"\t"<<psi<<"\t"<< distr->statisticalenergy(phi,psi)
				<<"\t"<< mixdistr->statisticalenergy(phi,psi)<<"\t"
				<<mixdistr->itplenergy(phi,psi,&dedp,&dedpsi)<<"\t";
		distr->randomphipsi(rng,&phi,&psi);
		std::cout << distr->statisticalenergy(phi,psi) <<"\t";
		double phir;
		double psir;
		distr->randompsi(rng,phi,&psir);
		std::cout <<distr->statisticalenergy(phi,psir)<<"\t";
		distr->randomphi(rng,&phir,psi);
		std::cout<<distr->statisticalenergy(phir,psi)<<std::endl;
	}*/
}


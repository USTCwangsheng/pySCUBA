/*
 * testssclassfier.cpp
 *
 *  Created on: 2017年12月10日
 *      Author: hyliu
 */
#include "pdbstatistics/nn_pepscorer.h"
#include "pdbstatistics/proteinblock.h"
#include "pdbstatistics/phipsidistr.h"
#include "dstl/randomengine.h"
using namespace NSPproteinrep;
using namespace NSPpdbstatistics;

int main(int argc,char **argv){
	std::vector<BackBoneSite> chain;
	readbackbonesites(std::string(argv[1]), chain);
	NN_SSClassifier &ssclassifier=NN_SSClassifier::getinstance();
	int length=7;

	for(auto it=chain.begin()+3;it !=chain.begin()+chain.size()-3;++it){
		if(!NSPproteinrep::fragstartsite(it-3,chain.end(),7)) continue;
		std::vector<double> pbtorsions=ProteinBlock::getpbtorsions(it,length);
		for(auto &t:pbtorsions) std::cout <<t<<" ";
		std::cout <<it->sscodechar()<<std::endl;
		std::vector<std::vector<double>> dp3dt;
		std::vector<double> p3=ssclassifier.probabilities(pbtorsions,&dp3dt);
/*		if(p3[0]>0.1 &&p3[0]<=0.5) {
			for(auto &t:pbtorsions) std::cout <<t<<" ";
			std::cout <<it->sscodechar()<<std::endl;
			for(int t=0;t<12;++t){
				std::cout <<" "<<dp3dt[0][t];
			}
			std::cout <<std::endl;
			std::vector<double> dp0dt(12,0.0);
			for(int m=0;m<pbtorsions.size();++m){
				pbtorsions[m] +=0.1;
				double p0p=ssclassifier.probabilities(pbtorsions,&dp3dt)[0];
				pbtorsions[m] -=0.2;
				double p0m=ssclassifier.probabilities(pbtorsions,&dp3dt)[0];
				dp0dt[m]=(p0p-p0m)/0.2;
				pbtorsions[m] +=0.1;
			}
			for(int t=0;t<12;++t){
				std::cout <<" "<<dp0dt[t];
			}
			std::cout <<std::endl<<std::endl;
		}*/
//		for(auto &t:pbtorsions) std::cout<<" "<<t;
//		std::cout <<std::endl;
		for(int s=0;s<3;++s) {
			std::cout <<"SS "<< s<<" prob: "<< p3[s];
			for(int t=0;t<12;++t){
				std::cout <<" "<<dp3dt[s][t];
			}
			std::cout <<std::endl;
		}
	}
/*	const PhiPsiDistr &distr=PhiPsiDistr::mixcoildistr();
	std::ofstream ofs;
	ofs.open("angcoil_r.dat");
	auto & rng=NSPdstl::RandomEngine<>::getinstance().realrng(0,1.0);
	for(int i=0;i<30000;++i) {
		std::vector<double> rtorsions;
		for(int l=0;l<length;++l){
			double phi,psi;
			distr.randomphipsi(rng,&phi,&psi);
			if(l==0) rtorsions.push_back(psi);
			else if(l==length-1) rtorsions.push_back(phi);
			else {rtorsions.push_back(phi);rtorsions.push_back(psi);}
		}
		for(auto &t:rtorsions) ofs <<t<<" ";
		ofs <<'C'<<std::endl;
	}*/
}

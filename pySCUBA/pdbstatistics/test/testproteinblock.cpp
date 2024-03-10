/*
 * testproteinblock.cpp
 *
 *  Created on: 2017年6月13日
 *      Author: hyliu
 */

#include "pdbstatistics/proteinblock.h"
#include "dstl/pca.h"
#include <iostream>
#include <fstream>
using namespace NSPpdbstatistics;
using namespace NSPproteinrep;
int main(int argc, char **argv) {
	std::vector<BackBoneSite> sites;
	readbackbonesites(std::string(argv[1]),sites);
	std::vector<char> pbtypes;
	std::vector<double> rmsdas;
	for (auto iter=sites.begin(); iter !=sites.end(); ++iter) {
		int posi=iter-sites.begin();
		char pb='x';
		double dev=0.0;
		if(posi >1 && fragstartsite(iter-2, sites.end(), 5,std::string(),false)){
			pb=ProteinBlock::pbtype(sites,posi,&dev);
		}
		pbtypes.push_back(pb);
		rmsdas.push_back(dev);
	}
	std::map<char,std::vector<std::vector<double>>> dihs;
	for(auto &pb:ProteinBlock::stdpbs)
		dihs.insert(std::make_pair(pb.first, std::vector<std::vector<double>>()));
	for(int i=0; i<pbtypes.size();++i) {
		if(pbtypes[i]=='x') {
			continue;
		}else {
			dihs[pbtypes[i]].push_back(ProteinBlock::extractdihvec(sites,i));
		}
	}
	for (auto & pbdih:dihs) {
		NSPdstl::PCA pca(8,true);
		pca.doPCA(pbdih.second);
		std::cout <<"PBTYPE: " <<pbdih.first<< '\t'<<pbdih.second.size()<<std::endl;
		for(int i=0; i<8;++i)
			std::cout <<'\t'<<ProteinBlock::stdpbs[pbdih.first][i];
		std::cout <<std::endl;
		const std::vector<double> & center=pca.getcenter();
		for(int i=0; i<8;++i)
				std::cout <<'\t'<<center[i];
		std::cout <<std::endl;
		std::vector<double> evec;
		std::vector<double> var;
		std::vector<double> acum;

		double var_tot=0;
		for(int i=0; i<8;++i) {
			evec.push_back(pca.geteigenvalue(i));
			if(i==0) acum.push_back(evec.back());
			else acum.push_back(acum.back()+evec.back());
			var.push_back(pca.getcovar(i,i));
			var_tot += var.back();
		}
		for(int i=0; i<8; ++i) {
			std::cout <<'\t' << evec[i];
		}
		std::cout <<std::endl;
		for(int i=0; i<8; ++i) {
			std::cout <<'\t' << acum[i]/var_tot;
		}
		std::cout <<std::endl;
		std::vector<double> cp1=pca.getcomponent(0);
		std::vector<double> cp2=pca.getcomponent(1);
		for(int i=0; i<8; ++i) {
			std::cout <<'\t' << cp1[i];
		}
		std::cout <<std::endl;
		for(int i=0; i<8; ++i) {
			std::cout <<'\t' << cp2[i];
		}
		std::cout <<std::endl;

		std::ofstream ofs;
		std::string filename="PB_x.dat";
		filename[3]=pbdih.first;
		ofs.open(filename.c_str());
		for(auto &dihvec:pbdih.second) {
			std::vector<double> pcacrd=pca.getPCAcrds(dihvec);
			for(int i=0; i<8; ++i)
				ofs <<' '<<pcacrd[i];
			ofs<<std::endl;
		}
		ofs.close();
	}
}

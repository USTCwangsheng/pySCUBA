/*
 * nn_pepscorer.cpp
 *
 *  Created on: 2017年12月10日
 *      Author: hyliu
 */
#include "pdbstatistics/nn_pepscorer.h"
using namespace NSPpdbstatistics;
const std::map<char,int> NN_SSClassifier::SSIDX{{'H',0},{'E',1},{'C',2}};
std::vector<double> NSPpdbstatistics::triangencoder(const std::vector<double> &torsions,
		std::vector<std::vector<double>> *dcdt){
	static std::vector<double> nang{1.,2.,4.,8.};
	static double degree=3.14159265/180.0;
	std::vector<double> code;
	if(dcdt){
		dcdt->assign(2*nang.size()*torsions.size(),std::vector<double>(torsions.size(),0.0));
	}
	for(int i=0;i<torsions.size();++i){
		double t=torsions[i]*degree;
		int cidx=2*i*nang.size();
		for(int n=0;n<nang.size();++n){
			code.push_back(cos(nang[n]*t));
			code.push_back(sin(nang[n]*t));
			if(dcdt){
				(*dcdt)[cidx+2*n][i]=-nang[n]*code[cidx+2*n+1]*degree;
				(*dcdt)[cidx+2*n+1][i]=nang[n]*code[cidx+2*n]*degree;
			}
		}
	}
	return code;
}
void NN_SSClassifier::setup(const std::string & filename){
//	std::string fn = NSPdataio::datapath() + filename;
	std::string fn=NSPdataio::datafilename(filename);
	std::ifstream ifs;
	ifs.open(fn.c_str());
//	char buffer[120];
//	ifs.getline(buffer,120);
//	encoder_intv_=std::stod(std::string(buffer));
	classifier_=NSPdstl::readl3softmaxclassifier(ifs);
	initialized_=true;
}
std::vector<double> NN_SSClassifier::probabilities(const std::vector<double> &pbtorsions,
		std::vector<std::vector<double>> *dp3dt){
	std::vector<std::vector<double>> dcdt;
	std::vector<double> tcodes=triangencoder(pbtorsions,&dcdt);
	std::vector<std::vector<double>> dp3dc;
	std::vector<double> p3=classifier_(tcodes,dp3dc);
	dp3dt->assign(3,std::vector<double>(pbtorsions.size(),0.0));
	int subtcodesize=tcodes.size()/pbtorsions.size();
	for(int i=0;i<3;++i){
		for(int j=0;j<subtcodesize;++j){
			for(int k=0;k<pbtorsions.size();++k){
				int jidx=k*subtcodesize+j;
				(*dp3dt)[i][k] +=dp3dc[i][jidx]*dcdt[jidx][k];
			}
		}
	}
	return p3;
}


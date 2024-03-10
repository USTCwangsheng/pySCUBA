/*
 * testligpsampler.cpp
 *
 *  Created on: 2020年10月24日
 *      Author: hyiu
 */
#include "sampling/ligpepsampler.h"
using namespace NSPintrct;
using namespace NSPsampling;
int main(int argc, char **argv){
	assert(argc>=2);
	std::string parfile(argv[1]);
	int nconfig=std::stoi(std::string(argv[2]));
	std::string jobname="auto";
	if(argc>3) jobname=std::string(argv[3]);
	IntrctPara ipara;
	SDRunPara sdpara;
	LigPSamplerPara lpspara;
	readpara_ligpepsampler(parfile,lpspara,ipara,sdpara,jobname);
	LigPepSampler sampler=mkligpepsampler(lpspara,ipara,sdpara);
	std::ofstream ligofs(lpspara.saveconfigfile);
	for(int i=0;i<nconfig;++i){
		sampler.buildandoptimize();
		ligofs <<i <<" ";
		sampler.saveconfig(ligofs);
	// for initial test
	//	std::string  pdb="ligconfig_"+std::to_string(i)+".pdb";
	//	std::ofstream opdb(pdb);
	//	sampler.imol()->writepdb(opdb);
	//	opdb.close();
	}
	ligofs.close();
}






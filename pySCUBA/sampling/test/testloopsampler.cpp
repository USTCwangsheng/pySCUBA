/*
 * loopsampling.cpp
 *
 *  Created on: 2019年12月23日
 *      Author: hyiu
 */

#include "sampling/loopsampler.h"

using namespace NSPintrct;
using namespace NSPsampling;
void printenergies(int lid, const std::array<double,IntrctBlck::ENESIZE> & energies,int loopidx=-1){
	std::string edescr="etot= ";
	if(loopidx>=0) edescr="eloop"+std::to_string(loopidx)+"= ";
	std::cout << "Loop " <<lid <<edescr;
	double etot=0.0;
	for(auto t:energies) etot+=t;
	std::cout<<etot << " e_components:  ";
	for(auto & t: energies) std::cout <<" "<<t;
	std::cout <<std::endl;
}
void writeloops(std::ostream &ofs,int lid, const std::array<double,IntrctBlck::ENESIZE> & energies,
		const std::vector<NSPgeometry::XYZ> & crds){
	ofs <<lid<<std::endl;
	for(auto & t:energies) ofs<<" "<<t;
	ofs<<std::endl;
	for(auto &c:crds) ofs <<c.toString()<<std::endl;
}
int main(int argc, char **argv){
	assert(argc>=2);
	std::string parfile(argv[1]);
	int nconfig=std::stoi(std::string(argv[2]));
	std::string jobname="auto";
	if(argc>3) jobname=std::string(argv[3]);
	IntrctPara ipara;
	SDRunPara sdpara;
	LpSamplerPara lpspara;
	readpara_lpsampler(parfile,lpspara,ipara,sdpara,jobname);
	LoopSampler sampler=mkloopsampler(lpspara,ipara,sdpara);
//	sampler.runMCCycles();
	double ene;
	std::vector<double> enes;
	std::array<double,IntrctBlck::ENESIZE> energies;
	std::vector<std::array<double,IntrctBlck::ENESIZE>> loopenergies;
	std::shared_ptr<std::vector<NSPgeometry::XYZ>> loopcrds;
	loopcrds=sampler.optimizeloops(ene,energies);
	int lidx=0;
	printenergies(lidx,energies);
    if(lpspara.nloops>1){
    	enes=sampler.loopenes(loopenergies);
    	for(int l=0;l<lpspara.nloops;++l){
    		printenergies(lidx,loopenergies[l],l);
    	}
    }
    std::ofstream nativeloops("startloops.dat");
    writeloops(nativeloops,lidx,energies,*loopcrds);
    nativeloops.close();
    std::ofstream outloops(lpspara.outloopfile);
	for(int i=0;i<nconfig;++i){
		lidx++;
		loopcrds=sampler.regenloops(ene,energies);
		printenergies(lidx,energies);
	    if(lpspara.nloops>1){
	    	enes=sampler.loopenes(loopenergies);
	      	for(int l=0;l<lpspara.nloops;++l){
	      		printenergies(lidx,loopenergies[l],l);
	      	}
	      }
		writeloops(outloops,lidx,energies,*loopcrds);
	}
}




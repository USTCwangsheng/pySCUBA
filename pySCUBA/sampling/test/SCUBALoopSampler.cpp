/*
 * SCUBALoopSampler.cpp
 *
 *  Created on: 2021年3月23日
 *      Author: yxchen
 */

#include "sampling/loopsampler.h"
#include "sampling/savedloops.h"
#include <random>
#include <cmath>
#include <fstream>
#include <string>
#include <iostream>
#include <stdlib.h>

using namespace NSPintrct;
using namespace NSPsampling;
using namespace NSPdataio;
void printenergies(int lid, const std::array<double,IntrctBlck::ENESIZE> & energies,int loopidx=-1){
	std::string edescr=" etot= ";
	if(loopidx>=0) edescr=" eloop"+std::to_string(loopidx)+"= ";
	std::cout << "Loop " <<lid <<edescr;
	double etot=0.0;
	for(auto t:energies)
		if (!isnan(t))
			etot+=t;
	std::cout<<etot << " e_components:  ";
	for(auto & t: energies)
		if (!isnan(t))
			std::cout <<" "<<t;
		else
			std::cout << " " << 0.0;
	std::cout <<std::endl;
}
void writeloops(std::string & nativeloops,int lid, const std::array<double,IntrctBlck::ENESIZE> & energies,
		const std::vector<NSPgeometry::XYZ> & crds){
	std::ofstream ofs(nativeloops, std::ofstream::app);
	ofs <<lid<<std::endl;
	for(auto & t:energies) ofs<<" "<<t;
	ofs<<std::endl;
	for(auto &c:crds) ofs <<c.toString()<<std::endl;
	ofs.close();
}

void baklspar() {
	std::ifstream f("bakls.par");
	if (f.good()) {
		system("mv bakls.par baklsold.par");
		system("mv ls.par bakls.par");
	}
	else {
		system("mv ls.par bakls.par");
	}
	f.close();
	std::string ls = "bakls.par";
	std::ifstream baklsfile(ls);
	std::ifstream newloopfile("newloop");
	std::string loopindx;
	getline(newloopfile, loopindx);
	std::string baklsline;
	std::string newls = "ls.par";
	std::ofstream newlsfile(newls);

	std::string::size_type idx;
	int n = 0;
	while (getline(baklsfile, baklsline)) {
		idx = baklsline.find("Loops");
		//std::cout << idx << std::endl;
		if (idx == -1) {
			newlsfile << baklsline;
			newlsfile << std::endl;
		}
		else {
				std::string loopline = "Loops = " + loopindx;
				newlsfile << loopline;
				newlsfile << std::endl;
			
		}
	}
	baklsfile.close();
	newloopfile.close();
	newlsfile.close();
}
void autoreadloop(const std::string &parfile) {
	baklspar();
	SavedLoops sl;
	ControlFile cf;
	cf.readfile(parfile);
	LpSamplerPara lpspara = makelpsamplerparam(std::string(), cf);
	sl.init(lpspara);
	sl.writermsdene("rmsdene.dat");
	std::ofstream ofs("nativequenched.pdb");
	sl.writepdb(ofs, sl.startconfig, true);
	int ntop = 100;
	if (lpspara.readloopnum > 0)
		ntop = lpspara.readloopnum;
	if (lpspara.readrmsd == "auto")
	{
		std::cout << "Will use " << sl.rmsdcut() * 10 << " Angstrom for RMSD_Cutoff" << std::endl;
		sl.savenonredundant("nonreduntloops.dat", "nr_rmsdene.dat", sl.rmsdcut());
		sl.writepdb_topn(ntop, sl.rmsdcut());
	}
	else if (lpspara.readrmsd == "min")
	{
		double nla = (double)sl.nloopatoms;
		double denominator = nla / sl.min_nloopatom;
		std::cout << "Will use " << sl.rmsdcut_de(denominator) * 10 << " Angstrom for RMSD_Cutoff" << std::endl;
		sl.savenonredundant("nonreduntloops.dat", "nr_rmsdene.dat", sl.rmsdcut_de(denominator));
		sl.writepdb_topn(ntop, sl.rmsdcut_de(denominator));
	}
	else
	{
		double rmsdcut = 0.0;
		//assert(rmsdcut = std::stod(lpspara.readrmsd));
		rmsdcut = atof(lpspara.readrmsd.c_str());
		std::cout << "Will use " << rmsdcut << " Angstrom for RMSD_Cutoff" << std::endl;
		sl.savenonredundant("nonreduntloops.dat", "nr_rmsdene.dat", rmsdcut / 10);
		sl.writepdb_topn(ntop, sl.rmsdcut_de(rmsdcut / 10));
	}
}

void searchloopandsample(LpSamplerPara& lpspara, IntrctPara& ipara, SDRunPara& sdpara) {
	LoopSampler sampler = mkloopsampler(lpspara, ipara, sdpara);
	sampler.readparams(lpspara);
	double ene, enestart, ene1, lastene, enegap;
	std::vector<double> enes;
	std::array<double, IntrctBlck::ENESIZE> energies;
	std::vector<std::array<double, IntrctBlck::ENESIZE>> loopenergies;
	std::shared_ptr<std::vector<NSPgeometry::XYZ>> loopcrds;
	//	loopcrds = sampler.loopatomcrds();
	loopcrds = sampler.optimizeloops(ene, energies); //pybind11
	//enestart = ene;
	int lidx = 0;
	printenergies(lidx, energies);
	if (lpspara.nloops > 1) {
		enes = sampler.loopenes(loopenergies);
		for (int l = 0;l < lpspara.nloops;++l) {
			printenergies(lidx, loopenergies[l], l);
		}
	}
	std::string nativeloops("startloops.dat");
	writeloops(nativeloops, lidx, energies, *loopcrds);
	//nativeloops.close();
	std::string outloops(lpspara.outloopfile);
	lastene = ene;
	sampler.imol()->getmolsystm().backupdata();
	// update sampler.looppools_ from startloop,in case cannot rebuild loop
	sampler.updatelooppools();

	std::cout << "Will Perform LoopSampling: " << lpspara.reconstructnum << " tries." << std::endl;
	for (int i = 0;i < lpspara.reconstructnum;++i) {
		//lidx++;
		//sampler.imol()->getmolsystm().backupdata();
		std::cout << "In reconstruct " << i << " round" << std::endl;
		loopcrds = sampler.regenuseloops(ene, energies);
		if (ene < lastene) {
			lidx++;
			std::cout << "new low ene :" << ene << " < lastene" << lastene << std::endl;
			lastene = ene;
			sampler.imol()->getmolsystm().backupdata();
			// update sampler.looppools_
			sampler.updatelooppools();
			printenergies(lidx, energies);
			if (lpspara.nloops > 1) {
				enes = sampler.loopenes(loopenergies);
				for (int l = 0;l < lpspara.nloops;++l) {
					printenergies(lidx, loopenergies[l], l);
				}
			}
			writeloops(outloops, lidx, energies, *loopcrds);
		}
		else {
			double p = ((double)rand()) / RAND_MAX;
			double expene = exp((lastene - ene) / 10);
			if (expene > p) {  //1 == rand() % 10
				lidx++;
				std::cout << "Accept a conformation with similar energy with a probability of " << expene << " > randomP: " << p << " , new ene :" << ene << "--lastene" << lastene << std::endl;
				lastene = ene;
				sampler.imol()->getmolsystm().backupdata();
				// update sampler.looppools_
				sampler.updatelooppools();
				printenergies(lidx, energies);
				if (lpspara.nloops > 1) {
					enes = sampler.loopenes(loopenergies);
					for (int l = 0;l < lpspara.nloops;++l) {
						printenergies(lidx, loopenergies[l], l);
					}
				}
				writeloops(outloops, lidx, energies, *loopcrds);
			}
			else {
				//When the energy does not meet the requirements, return to the previous conformation
				std::cout << "ene= " << ene << "> lastene= " << lastene << "  ,dropbackdata" << std::endl;
				sampler.imol()->getmolsystm().dropbackdata();
			}
		}
		std::cout << "lastene= " << lastene << std::endl;

	}
}
int main(int argc, char **argv){
	assert(argc>=2);
	std::string parfile(argv[1]);
	std::string jobname="auto";
	IntrctPara ipara;
	SDRunPara sdpara;
	LpSamplerPara lpspara;
	readpara_lpsampler(parfile,lpspara,ipara,sdpara,jobname);
	searchloopandsample(lpspara,ipara,sdpara);
	autoreadloop(parfile);
}

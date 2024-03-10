/*
 * SCUBALoopReader.cpp
 *
 *  Created on: 2021年3月23日
 *      Author: yxchen
 */

#include "sampling/savedloops.h"
using namespace NSPsampling;
using namespace NSPdataio;

LpSamplerPara loopreaderinit(SavedLoops& sl,std::string& parfile) {
	ControlFile cf;
	cf.readfile(parfile);
	LpSamplerPara lpspara = makelpsamplerparam(std::string(), cf);
	sl.init(lpspara);
	return lpspara;
}
void writestartpdb(SavedLoops& sl, std::string& startpdbconfig) {
	std::ofstream ofs(startpdbconfig);
	sl.writepdb(ofs, sl.startconfig, true);
}


void readloopout(SavedLoops& sl, LpSamplerPara& lpspara) {
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


int main(int argc, char **argv){
	SavedLoops sl;
	std::string parfile(argv[1]);
	LpSamplerPara lpspara = loopreaderinit(sl,parfile);
	sl.writermsdene("rmsdene.dat");
	std::string startpdbconfig="nativequenched.pdb";
	writestartpdb(sl, startpdbconfig);
	readloopout(sl, lpspara);
}

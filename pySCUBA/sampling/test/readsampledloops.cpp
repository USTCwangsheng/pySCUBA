/*
 * readsampledloops.cpp
 *
 *  Created on: 2019年12月25日
 *      Author: hyiu
 */

#include "sampling/savedloops.h"
using namespace NSPsampling;
using namespace NSPdataio;
int main(int argc, char **argv){
	SavedLoops sl;
	ControlFile cf;
	std::string parfile(argv[1]);
	cf.readfile(parfile);
	LpSamplerPara lpspara=makelpsamplerparam(std::string(),cf);
	sl.init(lpspara);
	sl.writermsdene("rmsdene.dat");
	std::ofstream ofs("nativequenched.pdb");
	sl.writepdb(ofs,sl.startconfig,true);
	sl.savenonredundant("nonreduntloops.dat","nr_rmsdene.dat",sl.rmsdcut());
	sl.writepdb_topn(100,sl.rmsdcut());
}



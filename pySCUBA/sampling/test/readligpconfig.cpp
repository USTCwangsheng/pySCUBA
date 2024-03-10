/*
 * readligpconfig.cpp
 *
 *  Created on: 2020年10月25日
 *      Author: hyiu
 */

#include "sampling/savedligpconfig.h"
using namespace NSPsampling;
using namespace NSPdataio;
int main(int argc, char **argv){
	SavedLigPConfig sl;
	ControlFile cf;
	std::string parfile(argv[1]);
	cf.readfile(parfile);
	LigPSamplerPara lpspara=makeligpsamplerparam(std::string(),cf);
	sl.init(lpspara);
	sl.writepdb_topn(100,0.3);
}




/*
 * testpocketmodeler.cpp
 *
 *  Created on: 2020年7月19日
 *      Author: hyiu
 */
#include "proteinrep/pocketmodeler.h"
using namespace NSPproteinrep;
int main(int argc,char **argv){
	PocketModeler modeler;
	modeler.interactivesetup(std::cout,std::cin);
	std::string file(argv[1]);
	std::ofstream outpdb0(file);
	modeler.writepdb(outpdb0);
	outpdb0.close();
	int nconf=10;
	for(int i=0;i<nconf;++i){
		bool success=modeler.genconformation(5000);
		if(!success) break;
		std::string filename=file+"_"+std::to_string(i)+".pdb";
		std::ofstream outpdb(filename);
		modeler.writepdb(outpdb);
	}
}




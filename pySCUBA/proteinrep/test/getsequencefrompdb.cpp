/*
 * getsequencefrompdb.cpp
 *
 *  Created on: 2020年8月24日
 *      Author: hyiu
 */
#include "proteinrep/pdbreader.h"
#include "proteinrep/aaconformer.h"
#include "proteinrep/aminoacidseq.h"
#include <iostream>
using namespace NSPproteinrep;

int main(int argc, char**argv){
	std::string filename(argv[1]);
	AAConformersInModel model;
	model.readpdbfile(filename);
	for(auto &c:model.conformers){
		std::vector<std::string> resnames;
		for(auto &r:c){
			resnames.push_back(r.residuename);
		}
		std::string codes=AminoAcidSeq::name2code(resnames);
		std::cout <<"<<chain-"<< c[0].chainid_or<<std::endl<<codes <<std::endl;
	}
	return 0;
}



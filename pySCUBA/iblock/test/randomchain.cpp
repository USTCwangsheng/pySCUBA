/*
 * randomchain.cpp
 *
 *  Created on: 2020年9月25日
 *      Author: hyiu
 */
#include "iblock/molmodeler.h"
#include "geometry/calculators.h"
#include "dstl/randomengine.h"
#include <dataio/splitstring.h>
void parsesseq(const std::string & sss, std::vector<char> &sseq){
	sseq.clear();
	std::vector<std::string> words;
	NSPutils::split(sss,words,":");
	for(auto &word : words){
        char c=word[0];
    	assert(c=='H' || c=='E' || c=='C');
    	std::vector<int> range=NSPdataio::integersInString(word.substr(1));
    	assert(range.size()==2);
		int len=NSPdstl::RandomEngine<>::getinstance().intrng(range[0],range[1])();
    	for(int i=0;i<len;++i){
    		sseq.push_back(c);
	}
	}
    return;
}
using namespace NSPintrct;
using namespace NSPproteinrep;
using namespace NSPgeometry;
int main(int argc, char**argv){
	const char *usage=
						R"""(usage:
        randomchain seed nconfiguration secstr 
        seed:  random number generator seed
        nconfiguration: number of configurations to generate
        secstr: specify type and length range of secondary structures along the chain, for example 
                        C1,30:H8,20:C5,30:H8,20:C5,30:H8,20:C5,30:H8,20:C5,30:H8,20:C1,30 
)""";
	if(argc !=4) {
		std::cout <<usage;
		exit(0);
	}
	int seed = std::stol(std::string(argv[1]));
	int nconfig = std::stol(std::string(argv[2]));
    std::string ssspec(argv[3]);
	NSPdstl::RandomEngine<>::getinstance().reseed(seed);
	MolModeler modeler;
	for (int nc = 0; nc < nconfig; ++nc) {
			IntrctMol imol;
			modeler.settarget(imol);
			std::vector<char> ssseq;
	        parsesseq(ssspec, ssseq);
			auto & rng = NSPdstl::RandomEngine<>::getinstance().realrng(0, 1);
			modeler.buildnewchain(XYZ(0,0,0), XYZ(1,0,0), ssseq);
			modeler.movechaincenter(XYZ(0,0,0),0);
			std::string filename = "config_" + std::to_string(nc) + ".pdb";
			std::ofstream ofs(filename);
			imol.writepdb(ofs);
			ofs.close();
			filename= "config_" + std::to_string(nc) +"_intendedss.txt";
			std::ofstream ofs1(filename);
			for(auto c:ssseq) ofs1<<c;
			ofs1<<std::endl;
			ofs1.close();
		}
}

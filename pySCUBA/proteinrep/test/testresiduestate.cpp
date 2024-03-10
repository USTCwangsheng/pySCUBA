/*
 * testresiduestates.cpp
 *
 *  Created on: 2020年12月29日
 *      Author: hyiu
 */

/*
 * This program shows how to use the ResidueState class
 * it reads a PDB file, calculate the secondary structure sequence, and output.
 */
#include "proteinrep/residuestate.h"
#include "proteinrep/aminoacidseq.h"
using namespace NSPproteinrep;
int main(int argc, char**argv){
	std::string pdbfile(argv[1]);
	std::ifstream is(pdbfile);
	auto rstates=residuestates(is);
	int chainidx=0;
	for(auto &c:rstates){
		std::string seq;
		std::string ssseq;
		std::cout<<"Chain " <<chainidx++<<":"<<std::endl;
		std::cout <<"seq\tresidue\tSS\tSAI\tPHI\tPSI"<<std::endl;
		int seqidx=1;
		for(auto &r:c){
			seq.push_back(AminoAcidSeq::name2code(r.sidechainstate.resiudetype));
			ssseq.push_back(r.backbonestate.SSState);
			std::cout<<seqidx++<<"\t"
							<< r.sidechainstate.resiudetype <<"\t"
					           <<r.backbonestate.SSState <<"\t"
							   << r.backbonestate.sai<<"\t"
							   <<r.backbonestate.phi<<"\t"
							   <<r.backbonestate.psi<<std::endl;
		}
		std::cout <<seq<<std::endl;
		std::cout<<ssseq<<std::endl;
	}
}

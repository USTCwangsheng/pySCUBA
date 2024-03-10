/*
 * preprocesspdb.cpp
 *
 *  Created on: 2020年10月23日
 *      Author: hyiu
 */
#include "iblock/intrctmol.h"
#include "proteinrep/aaconformer.h"
using namespace NSPintrct;
using namespace NSPproteinrep;

/**
 * this program preprocesses an input PDB file so that the output pdb can be used as input for
 * some other programs in  the SCUBA package.
 * 1. keeping only the protein atom records (MSE will replaced by MET)
 * 2. reassign chain IDs and renumber residues and atoms
 */
int main(int argc, char **argv) {
	std::string inpdbfile=std::string(argv[1]);
	std::string outpdbfile=std::string(argv[2]);
		AAConformersInModel model;
		std::shared_ptr<IntrctMol> imol;
	try {
		model.readpdbfile(inpdbfile);
	    imol=std::shared_ptr<IntrctMol>(new IntrctMol(model.conformers));
	}
	catch (...)
	{
		std::cout <<"Failed processing PDB file" <<std::endl;
		exit(1);
	}
	std::ofstream ofs(outpdbfile);
	ofs <<"PREPROCESSED PDB from " <<inpdbfile <<std::endl;
	imol->writepdb(ofs);
	ofs.close();
}




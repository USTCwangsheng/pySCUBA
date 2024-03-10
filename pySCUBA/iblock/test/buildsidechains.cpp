/*
 * addssbasedsidechains.cpp
 *
 *  Created on: 2020年10月2日
 *      Author: hyiu
 */
#include "proteinrep/aaconformer.h"
#include "proteinrep/aminoacidseq.h"
#include "iblock/molmodeler.h"
using namespace NSPproteinrep;
using namespace NSPintrct;
/*
 * build amino acid sidechains according to new user-supplied amino acid sequences or sidechain types mapped from
 * secondary structure types.
 * The secondary structure types can be automatically derived from the input structure or supplied by user
 */
int main(int argc, char **argv) {
	const char *usage=
						R"""(usage:
        buildsidechains inpdbfile outpdbfile  [aaseqfrom] [ssseqfrom] [helixtype strandtype coiltype] 
        aaseqfrom: "basedonss" or name of the file containing new one-letter aminoacid sequences
        ssseqfrom : "basedonstruct" or name of the file containined secondary structure sequences
        helixtype, strandtype, coiltype: mapping secondary structure type to residue type
)""";

	if(argc<3){
		std::cout<<usage<<std::endl;
		exit(0);
	}
	std::string inpdb(argv[1]);
	std::string outpdb(argv[2]);
	std::string aaseq { "basedonss" };
	if (argc >= 4) {
		aaseq = std::string(argv[4]); //"basedonss" or filename from which one letter code amino acid sequence will be read
	}
	std::string ssseq { "basedonstruct" };
	if (argc >= 5) {
		ssseq = std::string(argv[3]); // "fromstruct" or filename from which secondary struct sequence will be read
	}
	std::string helixtype { "LEU" };
	std::string strandtype { "VAL" };
	std::string coiltype { "GLY" };
	if (argc == 8) {
		helixtype = std::string(argv[5]);  //residue type for helix
		strandtype = std::string(argv[6]); //residue type for strand
		coiltype = std::string(argv[7]); //residue type for coil
	}
	AAConformersInModel model;
	model.readpdbfile(inpdb);
	std::shared_ptr<IntrctMol> imol = std::shared_ptr<IntrctMol>(
			new IntrctMol(model.conformers));
	std::vector<std::vector<std::string>> aatypes(imol->nchains());
	if (aaseq == "basedonss") {
		std::vector<std::string> secndstrs(imol->nchains(), std::string());
		if (ssseq == "basedonstruct") {
			secndstrs = imol->sscodestrings();
		} else {
			std::ifstream ifs(ssseq);
			for (int i = 0; i < imol->nchains(); ++i)
				ifs >> secndstrs[i];
			ifs.close();
		}
		for (int c = 0; c < imol->nchains(); ++c) {
			aatypes[c].resize(imol->nresidues(c));
			for (int r = 0; r < imol->nresidues(c); ++r) {
				if (secndstrs[c][r] == 'H')
					aatypes[c][r] = helixtype;
				else if (secndstrs[c][r] == 'E')
					aatypes[c][r] = strandtype;
				else
					aatypes[c][r] = coiltype;
			}
		}
	} else {
		std::ifstream ifs(aaseq);
		std::vector<std::string> seq1letter(imol->nchains());
		for (int c = 0; c < imol->nchains(); ++c){
			ifs >> seq1letter[c];
			assert(seq1letter[c].size() == imol->nresidues(c));
			aatypes[c] = AminoAcidSeq::code2name(seq1letter[c]);
		}
	}
	MolModeler modeler;
	modeler.checkclash(true);
	modeler.settarget(*imol);
	for (int c = 0; c < imol->nchains(); ++c) {
		for (int r = 0; r < imol->nresidues(c); ++r) {
			if (imol->residuename(NSPdstl::Idx2D(c, r)) != aatypes[c][r])
				modeler.changeresidue(NSPdstl::Idx2D(c, r), aatypes[c][r]);
		}
	}
	std::ofstream ofs(outpdb);
	imol->writepdb(ofs);
	ofs.close();
}

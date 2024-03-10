/*
 * SCUBAmakeLVG
 *
 *  Created on: 2019��12��18��
 *      Author: hyiu
 */
#include "iblock/molmodeler.h"
#include "iblock/intrctmol.h"
#include "iblock/intrctblck.h"
#include "designseq/StructureInfo.h"
#include "geometry/calculators.h"
#include "designseq/ProteinRep.h"
#include "dstl/randomengine.h"
#include "proteinrep/residuestate.h"
#include "proteinrep/aminoacidseq.h"

#include <dataio/splitstring.h>
#include <fstream>
#include <cassert>
#include <string>
#include <iostream> 

using namespace std;
using namespace NSPintrct;
using namespace NSPdesignseq;
using namespace NSPproteinrep;
using namespace NSPgeometry;



int main(int argc, char** argv) {
    if (argc != 3) {
        std::cout << "usage: SCUBAmakeLVG pdbfile outputfilename" << std::endl;
        exit(0);
    }

	std::string pdbfile(argv[1]);
	std::ifstream is(pdbfile);
	auto rstates = residuestates(is);
	int chainidx = 0;
	std::string ssseq;
	for (auto& c : rstates) {
		std::string seq;
		
		//std::cout << "Chain " << chainidx++ << ":" << std::endl;
		//std::cout << "seq\tresidue\tSS\tSAI\tPHI\tPSI" << std::endl;
		int seqidx = 1;
		for (auto& r : c) {
			seq.push_back(AminoAcidSeq::name2code(r.sidechainstate.resiudetype));
			ssseq.push_back(r.backbonestate.SSState);
			//std::cout << seqidx++ << "\t"
			//	<< r.sidechainstate.resiudetype << "\t"
			//	<< r.backbonestate.SSState << "\t"
			//	<< r.backbonestate.sai << "\t"
			//	<< r.backbonestate.phi << "\t"
			//	<< r.backbonestate.psi << std::endl;
		}
		std::cout << seq << std::endl;
		std::cout << ssseq << std::endl;
	}

int nres = ssseq.size();
MolModeler modeler;
//IntrctMol imol;
AAConformersInModel model;
model.readpdbfile(pdbfile);
std::shared_ptr<IntrctMol> imol = std::shared_ptr<IntrctMol>(new IntrctMol(model.conformers));
modeler.settarget(*imol);
modeler.checkclash(true);


for (int c = 0;c < nres;++c) {
    if (ssseq[c] == 'H') {
        modeler.changeresidue(MolModeler::SiteIdx{ 0,c }, "LEU"); //replace residue c (starting idx=0) of chain 0 by leu
    }
    else {
        if (ssseq[c] == 'E') {
            modeler.changeresidue(MolModeler::SiteIdx{ 0,c }, "VAL");
        }
        else {
            modeler.changeresidue(MolModeler::SiteIdx{ 0,c }, "GLY");
        }
    }
}
std::string out(argv[2]);
std::string filename = out + ".pdb";
std::ofstream ofs(filename);
imol->writepdb(ofs);
ofs.close();
std::string filename2 = out + "_LVGss.txt";
std::ofstream ofs1(filename2);
for (auto c : ssseq) ofs1 << c;
ofs1 << std::endl;
ofs1.close();


}

















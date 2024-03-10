/*
 * genrandomconfig.cpp
 *
 *  Created on: 2020年8月21日
 *      Author: hyiu
 */
#include "iblock/molmodeler.h"
#include "geometry/calculators.h"
#include "dstl/randomengine.h"
using namespace NSPintrct;
using namespace NSPproteinrep;
using namespace NSPgeometry;
/**
 * This program illustrates how to use the MolModeler class to generate
 * a cluster of peptide backbone fragments  (in strand or helix conformation) with random configurations
 * command line parameters are the following, in sequential orders
 * genrandomconfig  SEED NCONFIG NCHAINS LENMIN LENMAX SECSTRUC SPHERESIZE
 *SEED= seed for random number generator
 *NCONFIG= number of configurations to generate
 *NCHAINS=number of fragments to generate for each configuration
 *LENMIN=minimum length of fragments (for helix fragment, the minumm length will be 3*LENMIN)
 *LENMAX=maximum length of fragments (for helix fragment, the maximum length will be 3*LENMAX)
 *SECSTRUC=type of secondary structure, "E" for strand, "H" for helix, "R" for random choice between strand and helix
 *SPHERESIZE= N terminus of the fragments will be random points within a sphere of size SPHERESIZE
 */
#define PHELIX 0.3
int main(int argc, char **argv) {
	const char *usage =
			R"""(usage:
genrandomconfig  SEED NCONFIG NCHAINS LENMIN LENMAX SECSTRUC SPHERESIZE
SEED= seed for random number generator
NCONGIF=number of configurations to generate
NCHAINS=number of fragments to generate 
LENMIN=minimum length of fragments(for helix fragment, the minumm length will be 3*LENMIN)
LENMAX=maximum length of fragments(for helix fragment, the maximum length will be 3*LENMAX)
SECSTRUC=type of secondary structure,  "E" for strand,"H" for helix, "R" for random choice between strand and helix
SPHERESIZE= N terminus of the fragments will be random points within a sphere of size SPHERESIZE
)""";
	if (argc < 8) {
		std::cout << usage;
		exit(0);
	}
	int seed = std::stol(std::string(argv[1]));
	int nconfig = std::stol(std::string(argv[2]));
	int nchains = std::stol(std::string(argv[3]));
	int lenmin = std::stol(std::string(argv[4]));
	int lenmax = std::stol(std::string(argv[5]));
	char sstype = argv[6][0];
	assert(sstype == 'E' || sstype == 'H' || sstype == 'R');
	double radius = std::stod(std::string(argv[7]));
	NSPdstl::RandomEngine<>::getinstance().reseed(seed);
	MolModeler modeler;
	for (int nc = 0; nc < nconfig; ++nc) {
		IntrctMol imol;
		modeler.settarget(imol);
		for (int c = 0; c < nchains; c++) {
			auto & rng = NSPdstl::RandomEngine<>::getinstance().realrng(0, 1);
			char ss = sstype;
			if (sstype == 'R') {
				if (rng() < PHELIX)
					ss = 'H';
				else
					ss = 'E';
			}
			int lmin = lenmin;
			int lmax = lenmax;
			if (ss == 'H') {
				lmin *= 3;
				lmax *= 3;
			}
			auto &irng = NSPdstl::RandomEngine<>::getinstance().intrng(lmin,
					lmax);
			XYZ x1(rng, radius);
			double r1=sqrt(x1.squarednorm());
			XYZ x2(rng, radius);
			XYZ d = x2 - x1;
			d = 1.0 / (sqrt(d.squarednorm())) * d;
			int len = irng();
			std::vector<char> ssseq(len, ss);
			modeler.buildnewchain(x1, d, ssseq);
			modeler.movechaincenter(x1,c);
		}
		std::string filename = "config_" + std::to_string(nc) + ".pdb";
		std::ofstream ofs(filename);
		imol.writepdb(ofs);
		ofs.close();
	}
}


/*
 * testligandpep.cpp
 *
 *  Created on: 2020年10月22日
 *      Author: hyiu
 */
#include "sampling/ligandpep.h"
using namespace NSPgeometry;
using namespace NSPsampling;
using namespace NSPintrct;
using namespace NSPproteinrep;
int main(int argc, char **argv){
	   std::string pdbfilename(argv[1]);
		AAConformersInModel model;
		model.readpdbfile(pdbfilename);
		std::shared_ptr<IntrctMol> imol=std::shared_ptr<IntrctMol>(new IntrctMol(model.conformers));
		std::string selectionstr;
		for(int i=2;i<argc;++i) selectionstr=selectionstr +" "+ std::string(argv[i]);
		BlkSelector  sel(selectionstr);
		std::vector<NSPdstl::Idx2D> ifresidues;
		for(auto & entry:sel.selectedblks) {
			for(auto &res:entry.second){
				ifresidues.push_back(NSPdstl::Idx2D({entry.first,res}));
			}
		}
		LigandPep lp(*imol,ifresidues);
		int peplength=10;
		for(int m=0;m<20;++m){
			int ncid=imol->nchains();
			lp.addligandpep(imol.get(),peplength);
			lp.genligandpepcrd(imol.get(),ncid);
		}
		std::ofstream ofs("outligandpep.pdb");
		imol->writepdb(ofs);
}



/*
 * testmolmodeler.cpp
 *
 *  Created on: 2019年12月19日
 *      Author: hyiu
 */
#include "iblock/molmodeler.h"
#include "geometry/calculators.h"
#include "dstl/randomengine.h"
using namespace NSPintrct;
using namespace NSPproteinrep;
using namespace NSPgeometry;

/**
 * This program illustrates how to use the MolModeler class to change the structure of a protein,
 * including replacing side chains, rebuilding loops, split a polypeptide chain, joining
 * two chains into one with connecting loops, and so on
 *
 * Usage:
 * testmolmodeler inputpdbfile
 */
int main(int argc, char **argv){
#define BUILDCHAIN
#ifdef BUILDCHAIN
       int nchains=15;
       int lenmin=7;
       int lenmax=14;
       int seed=std::stol(std::string(argv[1]));
       NSPdstl::RandomEngine<>::getinstance().reseed(seed);
       IntrctMol imol;
       MolModeler modeler;
       modeler.settarget(imol);
       for (int c=0;c<nchains;c++){
    	   auto & rng=NSPdstl::RandomEngine<>::getinstance().realrng(0,1);
           auto &irng=NSPdstl::RandomEngine<>::getinstance().intrng(lenmin,lenmax);
    	   XYZ x1(rng,20.0);
    	   XYZ x2(rng,20.0);
    	   XYZ d=x2-x1;
    	   d=1.0/(sqrt(d.squarednorm()))*d;
    	   int len=irng();
    	   std::vector<char> ssseq(len,'E');
    	   modeler.buildnewchain(x1,d,ssseq);
       }
       std::ofstream ofs("randomfragments.pdb");
       imol.writepdb(ofs);
#endif
#ifdef CHANGECHAIN
	std::string pdbfilename(argv[1]);
	AAConformersInModel model;
	model.readpdbfile(pdbfilename);
	std::shared_ptr<IntrctMol> imol=std::shared_ptr<IntrctMol>(new IntrctMol(model.conformers));
	MolModeler modeler;
	modeler.checkclash(true);
	modeler.settarget(*imol);
	/**
	 * replace side chains
	 */
	modeler.changeresidue(MolModeler::SiteIdx{0,6},"GLU"); //replace residue 6 (starting idx=0) of chain 0 by glutamate
	modeler.changeresidue(MolModeler::SiteIdx{0,23},"PRO"); //replace residue 23 of chain 0 by proline

	/**
	 * rebuild loops
	 * setup confmaker for linkers between residue 5 and 12
	 */
	MolModeler::SiteIdx siteleft{0,15};
	MolModeler::SiteIdx siteright{0,21};
	auto confmaker=modeler.getloopconfmaker(siteleft,siteright);
	/**
	 * create a linker loop of 5 residues in length
	 */
	auto newloop=modeler.newloop(confmaker,std::vector<std::string>(5,"GLY"));
	if(!newloop.empty()){
		std::cout<<"Replacing a loop in target between  positions"
				<<siteleft.second<<" and "<<siteright.second<<std::endl;
		modeler.replaceloop(siteleft,siteright,newloop);
	}
	std::ofstream ofs("mutation.pdb");
	imol->writepdb(ofs);
	ofs.close();
	auto longerloop=modeler.newloop(confmaker,std::vector<std::string>(10,"GLY"));
		if(!newloop.empty()){
			std::cout<<"Intserting a longer loop in target between  positions"
					<<siteleft.second<<" and "<<siteright.second<<std::endl;
			auto newmol=modeler.mkimol_newloop(siteleft,siteright,longerloop);
			std::ofstream ofs("newmol.pdb");
			newmol.writepdb(ofs);
			ofs.close();
		}
		{
			std::cout<<"Split the molecule into two chains" <<std::endl;
			auto splitmol=modeler.splitchain(0,38,43);
			std::ofstream ofs("splitmol.pdb");
			splitmol.writepdb(ofs);
			ofs.close();

			MolModeler::SiteIdx left,right;
			modeler.settarget(splitmol);
			auto newmol1=modeler.mergechains(1,0,left,right);
			modeler.settarget(newmol1);
			ofs.open("mergemol.pdb");
			newmol1.writepdb(ofs);
			ofs.close();

			auto cfm=modeler.getloopconfmaker(left,right);
			auto cloop=modeler.newloop(cfm,std::vector<std::string>(12,"GLY"));
			if(!cloop.empty()){
				auto swapmol=modeler.mkimol_newloop(left,right,cloop);
				std::ofstream ofs("swapmol.pdb");
				swapmol.writepdb(ofs);
				ofs.close();
			}
		}
#endif
}

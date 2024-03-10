/*
 * ligandpep.h
 *
 *  Created on: 2020年10月17日
 *      Author: hyiu
 */

#ifndef LIGANDPEP_H_
#define LIGANDPEP_H_
#include "geometry/calculators.h"
#include "iblock/intrctmol.h"
namespace NSPsampling{

class LigandPep {
public:
	LigandPep(){;}
	LigandPep(const NSPintrct::IntrctMol & receptor,
	            const std::vector<NSPdstl::Idx2D> & interfaceresidues){
	            	setup(receptor,interfaceresidues);
	            }
	void setup(const NSPintrct::IntrctMol & receptor, const std::vector<NSPdstl::Idx2D> & interfaceresidues);
	std::vector<NSPgeometry::XYZ> genligandcrds(int ligandlength) const;
	NSPgeometry::XYZ randomligandcenter(
			NSPgeometry::XYZ *direction) const; //get a new center and direction for an initial ligand peptide
	void addligandpep(NSPintrct::IntrctMol *imol,int ligandlength,bool gencrd=false) const;
	void genligandpepcrd(NSPintrct::IntrctMol *imol,int ligandchainid) const;
private:
		std::vector<NSPgeometry::XYZ> interfacecacrds_;
		std::vector<std::pair<int,int>>  ifrespairs_;  //pairs of receptor residues at the interace, used for generate initial ligand peptide
		NSPgeometry::XYZ receptorcenter_;
};

/**
 * compare two configurations of a ligandpeptide by
 * sliding one along the other, and calculate the smallest RMSD of closest minaligment_ residues,
 * where minalignment_ is user provided or half the length of the peptide.
 * This can used to measure similarity between predicted configurations of a sequence-unspecialized
 *  peptide bound to a receptor
 */
class LigPConfigRMSD{
public:
	mutable int offsetmin{-1};
	LigPConfigRMSD(int peplength, int minalign=0):
		peplength_(peplength),minalignment_(minalign) {
		if(minalignment_==0) minalignment_=(peplength+1)/2;
	}
	double  operator()(const std::vector<NSPgeometry::XYZ> &cnfig1,
			const std::vector<NSPgeometry::XYZ> &cnfig2) const {
		double rmsdmin=1000000;
		offsetmin=-1;
		 for(int offset=0;offset<peplength_-minalignment_+1; ++offset){
			 std::vector<double> residuedevs;
			 for(int i=0;i<peplength_-offset;++i){
				  double diff2=0;
				  for(int m=0;m<4;++m){
					  NSPgeometry::XYZ dx=cnfig1.at(4*i+m)-cnfig2.at(4*i+m+4*offset);
					  diff2+=dx.squarednorm();
				  }
				  residuedevs.push_back(diff2);
			 }
			 std::sort(residuedevs.begin(),residuedevs.end());
			 double rmsd=0;
			 for(int i=0;i<minalignment_;++i){
				 rmsd +=residuedevs[i];
			 }
			 rmsd=sqrt(rmsd/(4.0*(double) minalignment_));
			 if(rmsd<rmsdmin) {
				 rmsdmin=rmsd;
				 offsetmin=offset;
			 }
		 }
		return rmsdmin;
	}
private:
	int minalignment_;
	int peplength_;
};
}


#endif /* LIGANDPEP_H_ */

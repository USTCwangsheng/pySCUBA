/*
 * rotamerdistributions.h
 *
 *  Created on: 2020年7月16日
 *      Author: hyiu
 */

#ifndef PROTEINREP_ROTAMERDISTRIBUTIONS_H_
#define PROTEINREP_ROTAMERDISTRIBUTIONS_H_
#include "designseq/ProteinRep.h"
#include "designseq/RotamerLib.h"
#include "dstl/randomengine.h"
#include <map>
namespace NSPproteinrep{
NSPdesignseq::PhipsiLib & getpplibinstance();
/**
 * phi-psi-independent probabilities of rotamers computed
 * according to phi-psi distributions and phi-psi-dependent rotamer energies
 */
struct RotamerDistribution{
	std::vector<NSPdesignseq::Rotamer *> rotamers;
	std::vector<double> rotamerene;
	std::vector<double> rotamerprob;
	/**
	 *setup the distribution of a residuetype, phi-psi-independent
	 *This is expected to be carried out at most once as needed.
	 */
	void setdistribution(const std::string & residuetype);

	/**
	 * setup the distribution of a residuetype, phi-psi-dependent.
	 */
	void setdistribution(const std::string &residuetype,double phi,double psi);
	static RotamerDistribution &getdistribution(const std::string &residuretype);
	/**
	 * draw a random rotamer according to the phi-psi-independent distribution;
	 * rrng() shoudl return a random ral number between 0 to 1.
	 */
	template<typename RNG>
	NSPdesignseq::Rotamer *randomrotamer(RNG &rng){
		while (true){
			int i=rng()*rotamers.size();
			if(i>rotamers.size()-1) i=rotamers.size()-1;
		    if(rng()<=rotamerprob[i]){
		    	return rotamers[i];
		    	}
		}
	}
};
inline NSPdesignseq::Rotamer *randomrotamer(const std::string &residuetype){
	RotamerDistribution & rd=RotamerDistribution::getdistribution(residuetype);
	auto &rng=NSPdstl::RandomEngine<>::getinstance().realrng(0,1);
	return rd.randomrotamer(rng);
}
inline NSPdesignseq::Rotamer *randomrotamer(const std::string &residuetype,double phi,double psi){
	RotamerDistribution rd;
	rd.setdistribution(residuetype,phi,psi);
	auto &rng=NSPdstl::RandomEngine<>::getinstance().realrng(0,1);
	return rd.randomrotamer(rng);
}
}




#endif /* PROTEINREP_ROTAMERDISTRIBUTIONS_H_ */

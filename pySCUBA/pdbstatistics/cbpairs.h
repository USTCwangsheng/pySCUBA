/*
 * cbpairs.h
 *
 *  Created on: 2018年1月24日
 *      Author: hyliu
 */

#ifndef PDBSTATISTICS_CBPAIRS_H_
#define PDBSTATISTICS_CBPAIRS_H_
#include "backbone/backbonesite.h"
namespace NSPpdbstatistics {

struct DisDistr{
	double step;
	int nsteps;
	std::vector<double> densities;
	DisDistr(double istep,int insteps):step(istep),nsteps(insteps){;}
	void estimate(std::vector<double> rsamples){
		std::vector<double> count(nsteps+1,0.0);
		double count_tot=0.0;
		for(auto r:rsamples){
			int idx=index(r);
			if(idx<nsteps){
				count[idx] +=1.0;
				count_tot+=1.0;
			}
		}
		double count_even=count_tot/(double) nsteps;
		for(int i=0;i<nsteps;i++){
			densities[i]=count[i]/count_even;
		}
		densities[nsteps]=1.0;
	}
	int index(double r) const{
		if (r<0) return 0;
		int idx=(int) (r/step);
		if (idx>nsteps) idx=nsteps;
		return idx;
	}
	double bincenter(int index) const {
		return ((double) index + 0.5)*step;
	}
	double densityintpl(double r) const {
		int i1=index(r-0.5*step);
		int i2=index(r+0.5*step);
		if(i1==i2) return densities[i1];
		double r1=bincenter(i1);
		double r2=bincenter(i2);
		double alpha=(r-r1)/r2-r1;
		return densities[i1]*alpha+densities[i2]*(1.0-alpha);
	}
/*	double eneintpl(double r) const {
		int i1=index(r-0.5*step);
		int i2=index(r+0.5*step);
		if(i1==i2) return -log(densities[i1]);
		double r1=bincenter(i1);
		double r2=bincenter(i2);
		double alpha=(r-r1)/r2-r1;
		return -log(densities[i1])*alpha-log(densities[i2])*(1.0-alpha);
	}*/
};

class PhiPsiRegion{
public:
	static std::set<int> resgionid;
	static int getregionid(double phi,double psi);
};

struct CBPair{
	CBPair(NSPproteinrep::BackBoneSite &bs1,NSPproteinrep::BackBoneSite &bs2);
	int phipsiresgion1_,phipsiregion2_;
	double rcb_{0.0};
};

CBPair randomcbpair(NSPproteinrep::BackBoneSite &bs1,NSPproteinrep::BackBoneSite &bs2);
}




#endif /* PDBSTATISTICS_CBPAIRS_H_ */

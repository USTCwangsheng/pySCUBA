/*
 * testgroupdistance.cpp
 *
 *  Created on: 2020年10月17日
 *      Author: hyiu
 */
/**
 * test groupdistance and its derivatives
 */
#include "sampling/groupdistance.h"
#include "dstl/randomengine.h"
#include <iostream>
using namespace NSPsampling;
using namespace NSPgeometry;

int main(int argc, char **argv){
	auto &rng=NSPdstl::RandomEngine<>::getinstance().realrng(0,1);
	NSPdstl::RandomEngine<>::getinstance().reseed(std::stoi(std::string(argv[1])));
	std::vector<XYZ> crd;
	for(int i=0;i<11;++i){
		crd.push_back(XYZ(rng,2.0));
	}
	GroupContact gc;
	auto &gd=gc.grpd;
//	gc.gdmin=std::stod(std::string(argv[2]));
//	gc.gdsmall=std::stod(std::string(argv[3]));
//	gc.gdoff=std::stod(std::string(argv[4]));
	std::vector<int> grp1{1,3,5},grp2{0,2,4,6,7,8,9,10};
	gd.grps[0]=grp1;
	gd.grps[1]=grp2;
	std::vector<NSPsampling::DvDxi> deriv;
	double d0=gc.theta(crd,&deriv);
	std::cout <<"Theta:: " <<d0 <<std::endl;
	std::vector<XYZ> dddr(crd.size(),{0,0,0});
	for(auto &d:deriv) dddr[d.first] =dddr[d.first]+d.second;
	for(auto i:grp1){
		for(int m=0;m<3; ++m){
			crd[i][m] +=0.0005;
			double dp=gc.theta(crd,&deriv);
			crd[i][m] -=0.001;
			double dm=gc.theta(crd,&deriv);
			crd[i][m]+=0.0005;
			std::cout <<i<<"  "<<m<<"  "<<(dp-dm)/0.001 <<"  .vs. " <<dddr[i][m]<<std::endl;
		}
	}
	for(auto i:grp2){
		for(int m=0;m<3; ++m){
			crd[i][m] +=0.0005;
			double dp=gc.theta(crd,&deriv);
			crd[i][m] -=0.001;
			double dm=gc.theta(crd,&deriv);
			crd[i][m]+=0.0005;
			std::cout <<i<<"  "<<m<<"  "<<(dp-dm)/0.001 <<"  .vs. " <<dddr[i][m]<<std::endl;
		}
	}
}


/*
 * groupdistance.cpp
 *
 *  Created on: 2020年10月17日
 *      Author: hyiu
 */
#include "sampling/groupdistance.h"

using namespace NSPsampling;
using namespace NSPgeometry;

double GroupDistance::distance(const std::vector<XYZ> & crd, std::vector<DvDxi> *dedxi) const {
		double r6sum=0.0;
		std::vector<NSPgeometry::XYZ> dedxi_1(grps[0].size(),{0,0,0});
		std::vector<NSPgeometry::XYZ> dedxi_2(grps[1].size(),{0,0,0});
		for(int i=0;i<grps[0].size();++i){
			for(int j=0;j<grps[1].size();++j){
				std::vector<NSPgeometry::XYZ> dedx;
				 double r=NSPgeometry::distance(crd[grps[0][i]],crd[grps[1][j]],&dedx);
				 r6sum+=1/pow(r,6);
				 NSPgeometry::XYZ dr6dri=-6/pow(r,7)*dedx[0];
				 dedxi_1[i]=dedxi_1[i]+dr6dri;
				 dedxi_2[j]=dedxi_2[j]-dr6dri;
			}
		}
//		double npairs=(double) (grps[0].size()*grps[1].size());
//		double r6sumav=r6sum/npairs;
//		double d=1.0/pow(r6sumav,1./6.0);
		double d=1.0/pow(r6sum,1.0/6.0);
//		double drdr6sumav=-(1./6.0)/pow(r6sumav,7.0/6.0)/npairs;
		double drdr6sum=-(1.0/6.0)/pow(r6sum,7.0/6.0);
		dedxi->clear();
		for(int i=0;i<grps[0].size();++i){
//			dedxi->push_back(std::make_pair(grps[0][i],drdr6sumav*dedxi_1[i]));
			dedxi->push_back(std::make_pair(grps[0][i],drdr6sum*dedxi_1[i]));
		}
		for(auto i=0;i<grps[1].size();++i){
	//		dedxi->push_back(std::make_pair(grps[1][i],drdr6sumav*dedxi_2[i]));
			dedxi->push_back(std::make_pair(grps[1][i],drdr6sum*dedxi_2[i]));
		}
//		std::cout <<"Group distance: " <<d <<std::endl;
		return d;
	}
double GroupContact::theta(const std::vector<NSPgeometry::XYZ> &crd,std::vector<DvDxi> *dthetadxi) const {
		std::vector<DvDxi> dgddxi;
		double dthetadgd;
		double res=theta(grpd.distance(crd,&dgddxi),&dthetadgd);
		if(dthetadgd==0.0) return res;
		for(auto &d:dgddxi){
			dthetadxi->push_back(std::make_pair(d.first,dthetadgd*d.second));
		}
		return res;
	}
	double GroupContact::theta(double gd,double *deriv) const {
        *deriv=0.0;
//		std::cout << "gdmin=" << gdmin << ", gdsmall=" << gdsmall << ", gdoff=" << gdoff << std::endl;
		if(gd<gdmin) return 1.0;
		if(gd>gdoff) return 0.0;
		double res;
		if(gd<gdsmall){
			 res=1.0-(1.0-thetasmall)*(gd-gdmin)/(gdsmall-gdmin);
			 *deriv=-(1.0-thetasmall)/(gdsmall-gdmin);
		} else {
			res= thetasmall-thetasmall*(gd-gdsmall)/(gdoff-gdsmall);
			*deriv=-thetasmall/(gdoff-gdsmall);
		}
		return res;
	}


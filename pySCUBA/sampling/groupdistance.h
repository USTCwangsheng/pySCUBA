/*
 * groupdistance.h
 *
 *  Created on: 2020年10月16日
 *      Author: hyiu
 */

#ifndef GROUPDISTANCE_H_
#define GROUPDISTANCE_H_
#include "geometry/calculators.h"
#include <array>
#include <vector>
namespace NSPsampling{
typedef std::pair<int, NSPgeometry::XYZ> DvDxi;
/**
 * approximate minimum inter-atomic distances between two group of atoms
 * distance=(average (rij^(-6))^(-1/6)), for all i-j atom pairs between grps 1 and 2.
 */
struct GroupDistance{
	std::array<std::vector<int>,2> grps; //atoms of each group
	double distance(const std::vector<NSPgeometry::XYZ> &crd, std::vector<DvDxi> *drdxi) const;
};

/**
 * determine if two groups of atoms are in contact.
 * theta returns 1 if group distance below gdmin, 0 if group distance above gdmax
 * and between 0 to 1 for group distance between gdmin and gdmax
 */
struct GroupContact{
	GroupDistance grpd;
	double gdmin{0.9},gdsmall{1.0},gdoff{3.5}; //theta=1 for gd below gdmin,
	                                                            //linearly from 1 to theasmall from gdmin to gdsmall,
	                                                           // and from thetasmall to 0 from gdsmall to gdoff
	double thetasmall{0.1};
	double theta(const std::vector<NSPgeometry::XYZ> &crd,std::vector<DvDxi> *dthetadxi) const;
	double theta(double gd,double *deriv) const;
};
}




#endif /* GROUPDISTANCE_H_ */

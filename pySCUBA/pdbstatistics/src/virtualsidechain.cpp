/*
 * virtualsidechain.cpp
 *
 *  Created on: 2018年2月3日
 *      Author: hyliu
 */



#include "pdbstatistics/virtualsidechain.h"

using namespace NSPgeometry;
using namespace NSPpdbstatistics;
using namespace NSPproteinrep;
std::map<std::string,NSPgeometry::XYZ> VirtualSideChain::localvcrd{
{"VB", {-0.484,-0.738,1.249}},
{"VG1", {0.038,-2.175, 1.228}},
{"VG2", { -2.014, -0.752,  1.273}},
{"VD1", {  -0.447,-2.913,   2.477}},
{"VD2",{  -2.498,  -1.490,   2.522}}
};


VirtualSideChain::VirtualSideChain(const XYZ &cacrd,const XYZ &ncrd,const XYZ &ccrd, double lengthscale):
	a3localframe_(cacrd,ncrd,ccrd){
	for(auto &lvcd:localvcrd){
		vcrds_.insert(std::make_pair(lvcd.first,a3localframe_.local2globalcrd(lengthscale*lvcd.second)));
	}
/*	static const double deg=3.14159265/180.0;
	NSPgeometry::XYZ rvb=InternaltoXYZ(cacrd,ccrd,ncrd,1.53,109.5*deg,120.0*deg);
	NSPgeometry::XYZ rvg1=InternaltoXYZ(rvb,cacrd,ncrd,1.53,109.5*deg,-60*deg);
	NSPgeometry::XYZ rvg2=InternaltoXYZ(rvb,cacrd,ncrd,1.53,109.5*deg,180*deg);
	NSPgeometry::XYZ rvd1=InternaltoXYZ(rvg1,rvb,cacrd,1.53,109.5*deg,180.0*deg);
	NSPgeometry::XYZ rvd2=InternaltoXYZ(rvg2,rvb,cacrd,1.53,109.5*deg,180.0*deg);
	vcrds_.insert(std::make_pair(std::string("VB"),rvb));
	vcrds_.insert(std::make_pair(std::string("VG1"),rvg1));
	vcrds_.insert(std::make_pair(std::string("VG2"),rvg2));
	vcrds_.insert(std::make_pair(std::string("VD1"),rvd1));
	vcrds_.insert(std::make_pair(std::string("VD2"),rvd2));*/
}
std::vector<XYZ> VirtualSideChain::redistributederiv(const std::string &atomname,
		const NSPgeometry::XYZ & deriv_vc) const {
		return a3localframe_.distributedvdx(vcrds_.at(atomname),deriv_vc);
}

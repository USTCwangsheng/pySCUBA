/*
 * virtualsidechain.h
 *
 *  Created on: 2018年2月3日
 *      Author: hyliu
 */

#ifndef PDBSTATISTICS_VIRTUALSIDECHAIN_H_
#define PDBSTATISTICS_VIRTUALSIDECHAIN_H_

#include "backbone/backbonesite.h"
#include "geometry/localframe.h"
#include <map>
namespace NSPpdbstatistics{

class VirtualSideChain{
public:
	static std::map<std::string,NSPgeometry::XYZ> localvcrd;
	VirtualSideChain(){;}
	VirtualSideChain(const NSPgeometry::XYZ & cacrd,
			const NSPgeometry::XYZ &ncrd, const NSPgeometry::XYZ &ccrd,double lengthscale=1.0);
	const NSPgeometry::XYZ & getcrd(const std::string &atomname )const{
		return vcrds_.at(atomname);
	}
	NSPgeometry::XYZ getlocalcrd(const std::string &atomname )const{
		return a3localframe_.global2localcrd(vcrds_.at(atomname));
	}
	std::vector<NSPgeometry::XYZ> redistributederiv(const std::string &atomname,
			const NSPgeometry::XYZ & deriv_vc) const ;

private:
	std::map<std::string,NSPgeometry::XYZ> vcrds_;
	NSPgeometry::A3LocalFrame a3localframe_;
};

}



#endif /* PDBSTATISTICS_VIRTUALSIDECHAIN_H_ */

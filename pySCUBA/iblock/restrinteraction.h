/*
 * restraintrction.h
 *
 *  Created on: 2020年7月22日
 *      Author: hyiu
 */

#ifndef IBLOCK_RESTRAINTRCTION_H_
#define IBLOCK_RESTRAINTRCTION_H_
#include "iblock/intrctmol.h"
#include "iblock/blcktopo.h"
namespace NSPintrct{
class RestrainFunc{
public:
	double operator()(double x, double &dedx) const{
		double x0;
		if(x<xmin_) x0=xmin_;
		else if(x>xmax_) x0=xmax_;
		else {
			dedx=0;
			return x;
		}
		double diff=x-x0;
		dedx=k_*diff;
		return 0.5*dedx*diff;
	}
	RestrainFunc(double xmin=-1.e6, double xmax=1.e6,double k=0.0):
		xmin_(xmin),xmax_(xmax),k_(k){;}
	//todo implement the following methods
	double distanceforce(const NSPgeometry::XYZ & x1,
			const NSPgeometry::XYZ &x2, std::vector<NSPgeometry::XYZ> *dedx);
	double angleforce(std::vector<NSPgeometry::XYZ> &crds,
			std::vector<NSPgeometry::XYZ> *dedx);
	double torsionforce(std::vector<NSPgeometry::XYZ> &crds,
				std::vector<NSPgeometry::XYZ> *dedx);
private:
	double xmin_;
	double xmax_;
	double k_;
};
class RestrInteraction{
public:
	virtual double force(IntrctMol &molsys,std::vector<NSPgeometry::XYZ>  &dedx)=0;
	virtual ~RestrInteraction(){;}
};

class PocketGeomRestr:public RestrInteraction{
public:
	//todo implement the force method
	virtual double force(IntrctMol &molsys,std::vector<NSPgeometry::XYZ>  &dedx);
	void setjointatom(int i, const IntrctBlck &blck, const std::string &atomname){
		int aid=blck.topo->atomindex(atomname);
		assert(aid>=0);
	    setjointatom(i,blck.getidx2d().first,blck.getidx2d().second,aid);
	}
	void setjointatom(int i, int chainid, int resid,int aid){
		jointatoms_[i][0]=chainid;
		jointatoms_[i][1]=resid;
		jointatoms_[i][2]=aid;
	}
	void restraintorsion0123(double xmin,double xmax,double k ){
		restrainfuncs_[0]=RestrainFunc(xmin*3.14159265/180.0,
				xmax*3.14159265/180.0,k);
	}
	void restrainangle123(double xmin,double xmax,double k ){
			restrainfuncs_[1]=RestrainFunc(xmin*3.14159265/180.0,
					xmax*3.14159265/180.0,k);
		}
	void restraindistance23(double xmin,double xmax,double k ){
				restrainfuncs_[2]=RestrainFunc(xmin,
						xmax,k);
			}
	void restraintorsion1234(double xmin,double xmax,double k ){
		restrainfuncs_[3]=RestrainFunc(xmin*3.14159265/180.0,
				xmax*3.14159265/180.0,k);
	}
	void restrainangle234(double xmin,double xmax,double k ){
				restrainfuncs_[4]=RestrainFunc(xmin*3.14159265/180.0,
						xmax*3.14159265/180.0,k);
			}
	void restraintorsion2345(double xmin,double xmax,double k ){
				restrainfuncs_[5]=RestrainFunc(xmin*3.14159265/180.0,
						xmax*3.14159265/180.0,k);
			}
private:
	std::array<AtomIdx3D, 6> jointatoms_;
	std::array<RestrainFunc,6> restrainfuncs_;
};
}


#endif /* IBLOCK_RESTRAINTRCTION_H_ */

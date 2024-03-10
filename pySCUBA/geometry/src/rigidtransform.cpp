/*
 * rigidtransform.cpp
 *
 *  Created on: 2016年11月8日
 *      Author: hyliu
 */

#include "geometry/rigidtransform.h"
#include "dstl/randomengine.h"
#include "geometry/quatfit.h"

using namespace NSPgeometry;

RigidTransform NSPgeometry::applyRotation(const Rotation &r, const RigidTransform &rt){
	RigidTransform rtnew;
	typename Rotation::constmatrixpointer m=rt.rotation().matrix();
	typename Rotation::constmatrixpointer rm=r.matrix();
	Rotation & rnew=rtnew.rotation();
	XYZ & trans=rtnew.translation();
	typename Rotation::matrixpointer mnew=rnew.matrix();
	for(int i=0;i<3;++i) {
		for(int j=0; j<3;++j)
			mnew[i][j]= rm[i][0]*m[0][j]+rm[i][1]*m[1][j]+rm[i][2]*m[2][j];
	}
	rnew.center()= rt.rotation().center();
	trans= rt.translation()+rnew.center();
	r.apply(&trans);
	trans= trans -rnew.center();
	return rtnew;
}

void RigidTransform::apply(XYZ *p) const{
	rot_.apply(p);
	*p = *p + trans_;
}
XYZ RigidTransform::applytoCopy(const XYZ &p) const {
	XYZ res=p;
	apply(&res);
	return res;
}

RigidTransform::RigidTransform(const XYZ & ri, const XYZ &rj, const XYZ &rk, const XYZ &rl,
		double dkl, double tjkl, double pijkl){
        XYZ rkl=rl-rk;
        XYZ rjk=rk-rj;
        XYZ rjkl=rk+cross(rjk,rkl);
        double d0kl=sqrt(rkl.squarednorm());
        Rotation rott(rj,rk,rjkl,rl,-tjkl);
        Rotation rotp(ri,rj,rk,rl,pijkl);
        *this=rotp*rott;
        XYZ t=(dkl-d0kl)/d0kl*(this->applytoCopy(rl)-rk);
        this->translation()=this->translation()+t;
}
RigidTransform::RigidTransform(const std::vector<XYZ> &jointpts, const std::vector<double> &intcrds){
	RigidTransform t0(jointpts[0],jointpts[1],jointpts[2],jointpts[3],intcrds[2],intcrds[1],intcrds[0]);
	RigidTransform t1(Rotation(jointpts[1],jointpts[2],t0.applytoCopy(jointpts[3]),t0.applytoCopy(jointpts[4]),intcrds[3]));
	RigidTransform t10=t1*t0;
	RigidTransform rt2(t10.applytoCopy(jointpts[5]),t10.applytoCopy(jointpts[4]),
			t10.applytoCopy(jointpts[3]),jointpts[2],intcrds[2],intcrds[4],intcrds[5]);
	*this=rt2.getreverse()*t10;
}

RigidTransform NSPgeometry::operator *(const RigidTransform &  a, const RigidTransform &b){
			RigidTransform rt;
			auto mt=rt.rotation().matrix();
			RigidTransform a0=a.axistoorigin();
			RigidTransform b0=b.axistoorigin();
			auto ma=a0.rotation().matrix();
			auto mb=b0.rotation().matrix();
			for(int i=0;i<3;++i){
				for(int j=0;j<3;++j){
					mt[i][j]=ma[i][0]*mb[0][j]+ma[i][1]*mb[1][j]+ma[i][2]*mb[2][j];
				}
			}
			rt.translation()=a0.rotation().applytoCopy(b0.translation())+a0.translation();
			return rt;
}
RigidTransform NSPgeometry::randomrigidtransform(double maxrotate, double maxtranslate){
	static auto &reng=NSPdstl::RandomEngine<>::getinstance();
	double angle=maxrotate*(reng.realrng(0.0,1.0)()-0.5);
	XYZ axis(reng.realrng(0.0,1.0),1.0);
	return RigidTransform(QuaternionCrd(axis,angle),XYZ(0.0,0.0,0.0),
				XYZ(reng.realrng(0.0,1.0),maxtranslate));
}

RigidTransform NSPgeometry::superpose(const std::vector<XYZ> &crda, const std::vector<XYZ> &crdb,
		const std::vector<std::pair<int,int>> & alignedpositions,double *rmsd2){
		std::vector<XYZ> refcrd;
		std::vector<XYZ> crd;
		for(auto & ap:alignedpositions){
			refcrd.push_back(crda[ap.first]);
			crd.push_back(crdb[ap.second]);
		}
		QuatFit qf;
		double dev2=qf.setup(refcrd,crd);
		if(rmsd2) *rmsd2=dev2;
		return qf.getRigidTransform();
}


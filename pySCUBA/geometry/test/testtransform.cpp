/*
 * testrotation.cpp
 *
 *  Created on: 2016年9月29日
 *      Author: hyliu
 */


#include "geometry/rigidtransform.h"
#include "geometry/quatfit.h"
using namespace NSPgeometry;
int main(int argc,char **argv){
/*	Rotation r1(QuaternionCrd(XYZ(1,2,3),30.0),XYZ());
	RigidTransform rt(QuaternionCrd(XYZ(1,2,3),30.0),XYZ(),XYZ(3,-1,4));
        auto q=QuatFit();
	XYZ pa(3.2,4.3,1.5);
	XYZ pbb=r1.applytoCopy(pa);
	XYZ pb=rt.applytoCopy(pa);
	std::cout <<pb.toString() <<std::endl;
	std::cout <<pbb.toString() <<std::endl;
	Rotation r2(QuaternionCrd(XYZ(1,0,-3),60.0),XYZ(1,3,5));
	RigidTransform rtnew=applyRotation(r2,rt);
	XYZ pc=r2.applytoCopy(pb);
	XYZ pcc=rtnew.applytoCopy(pa);
	std::cout <<pc.toString() <<std::endl;
	std::cout <<pcc.toString() <<std::endl;
*/
/** six pts, the first three are ifx,jfx, kfx, the second three are kmv,jmv,imv
 *
 */
	std::vector<XYZ>  pts{
		{1,2,3},{1,2,4},{1,4,4},
		{1,4,5},{2,4,5},{3,5,5}
	};
	/**
	 * six target internal coordinate values. angles are in radian
	 */
	std::vector<double> intcrds{
		-2.0,  2.8,1.6,-3.1,2.5,1.9
	};
	RigidTransform rt(pts,intcrds);
	for(int i=3;i<6;++i){
		rt.apply(&(pts[i]));
	}

	/**
	 * compare the actual and target internal coordinate values
	 */
	std::cout <<"torsion (-2.0) "<< torsion(pts[0],pts[1],pts[2],pts[3])<<std::endl;
	std::cout <<"angle(2.8)  "<< angle(pts[1],pts[2],pts[3])<<std::endl;
	std::cout <<"distance (1.6)  " << distance(pts[2],pts[3])<<std::endl;
	std::cout <<"torsion (-3.1) "<< torsion(pts[1],pts[2],pts[3],pts[4])<<std::endl;
	std::cout <<"angle(2.5)  "<< angle(pts[2],pts[3],pts[4])<<std::endl;
	std::cout <<"torsion (1.9) "<< torsion(pts[2],pts[3],pts[4],pts[5])<<std::endl;
/*	XYZ p1(3.0,-5.4,8.3);
	XYZ p2(4.0,-7.4,9.3);
	XYZ f1(-2,3,0.5);
	XYZ f2(0,4.1,-0.5);
	Rotation R=rotationaligntwovectors(p1,p2-p1,f2-f1);
	XYZ newf2=p1+f2-f1;
	R.apply(&newf2);
	std::cout <<"p2: " <<(p2-p1).toString() <<std::endl;
	std::cout <<"newf2: " <<(newf2-p1).toString() <<std::endl;

	XYZ r1(0.2,0.3,0.5);
	double theta1=sqrt(r1.squarednorm())/3.14159265*180.0;
	QuaternionCrd q1(r1,theta1);
	XYZ r2(0.1,-0.2,0.4);
	double theta2=sqrt(r2.squarednorm())/3.14159265*180.0;
	QuaternionCrd q2(r2,theta2);
	QuaternionCrd q3=q2*q1;
	std::cout <<q3.diff(q2) <<"   " <<theta1<<std::endl;
	Rotation m;
	m.init(q2,XYZ(0,0,0));
	std::cout <<r1.toString() <<std::endl;
	m.apply(&r1);
	std::cout <<r1.toString() <<std::endl;
	XYZ rt=r1+r2;
	double theta4=sqrt(rt.squarednorm())/3.14159265*180.0;
	QuaternionCrd q4(rt,theta4);
	std::cout <<q3.diff(q4)<<std::endl;
	*/
}


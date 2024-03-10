/*
 * loopsolver.cpp
 *
 *  Created on: 2016年12月7日
 *      Author: hyliu
 */
#include "loopclosure/loopsolver.h"
#include "loopclosure/loopclosure.h"
#include "geometry/calculators.h"
#include "geometry/quatfit.h"
#include <iostream>
using namespace NSPloopclosure;
using namespace NSPgeometry;


unsigned int LoopSolver::getSolutions(std::vector<XYZ> & atomcrds,
		std::vector<LoopSolution> & solutions){
	std::vector<XYZ> fixedcrds;
	fixedcrds.push_back(atomcrds[N1]);
	fixedcrds.push_back(atomcrds[CA1]);
	fixedcrds.push_back(atomcrds[CA3]);
	fixedcrds.push_back(atomcrds[C3]);
	return getSolutions(atomcrds,fixedcrds,solutions);
}
unsigned int LoopSolver::getSolutions(std::vector<XYZ> & atomcrds, std::vector<XYZ> &fixedcrds,
		std::vector<LoopSolution> & solutions) {
	double b_len[6], b_ang[7], t_ang[2];
	for (unsigned int i = 1; i < 7; ++i)
		b_len[i-1] = sqrt((atomcrds[i + 1] - atomcrds[i]).squarednorm());
	for (unsigned int i = 0; i < 7; ++i) {
		b_ang[i] = angle(atomcrds[i], atomcrds[i + 1], atomcrds[i + 2]);
		if(b_ang[i] < 0.000001) return 0;
	}
	t_ang[0] = torsion(atomcrds[C1], atomcrds[CA1], atomcrds[CA2],
			atomcrds[N2]);
	t_ang[1] = torsion(atomcrds[C2], atomcrds[CA2], atomcrds[CA3],
			atomcrds[N3]);
	initialize_loop_closure(b_len, b_ang, t_ang);
	double r_n1[3], r_a1[3], r_a3[3], r_c3[3];
	r_n1[0] = fixedcrds[0].x_;
	r_n1[1] = fixedcrds[0].y_;
	r_n1[2] = fixedcrds[0].z_;
	r_a1[0] = fixedcrds[1].x_;
	r_a1[1] = fixedcrds[1].y_;
	r_a1[2] = fixedcrds[1].z_;
	r_a3[0] = fixedcrds[2].x_;
	r_a3[1] = fixedcrds[2].y_;
	r_a3[2] = fixedcrds[2].z_;
	r_c3[0] = fixedcrds[3].x_;
	r_c3[1] = fixedcrds[3].y_;
	r_c3[2] = fixedcrds[3].z_;
	const int max_soln=16;
	int n_soln;
//
	static int ncount{0};
//	std::cout <<"Ncalls "<<ncount++ <<std::endl;
	double r_soln_n[max_soln][3][3], r_soln_a[max_soln][3][3], r_soln_c[max_soln][3][3];
	solve_3pep_poly(r_n1, r_a1, r_a3,
			r_c3, r_soln_n, r_soln_a,
			r_soln_c, &n_soln);
	if(n_soln<=0) return n_soln;
	std::vector<XYZ> block1;
	std::vector<XYZ> block2;
	block1.push_back(atomcrds[CA1]);
	block1.push_back(atomcrds[C1]);
	block1.push_back(atomcrds[N2]);
	block1.push_back(atomcrds[CA2]);
	block2.push_back(atomcrds[CA2]);
	block2.push_back(atomcrds[C2]);
	block2.push_back(atomcrds[N3]);
	block2.push_back(atomcrds[CA3]);
	for(unsigned int k=0; k<n_soln;++k) {
		std::vector<XYZ> refblock1;
		std::vector<XYZ> refblock2;
		refblock1.push_back(XYZ(r_soln_a[k][0][0],r_soln_a[k][0][1],r_soln_a[k][0][2]));
		refblock1.push_back(XYZ(r_soln_c[k][0][0],r_soln_c[k][0][1],r_soln_c[k][0][2]));
		refblock1.push_back(XYZ(r_soln_n[k][1][0],r_soln_n[k][1][1],r_soln_n[k][1][2]));
		refblock1.push_back(XYZ(r_soln_a[k][1][0],r_soln_a[k][1][1],r_soln_a[k][1][2]));
		refblock2.push_back(XYZ(r_soln_a[k][1][0],r_soln_a[k][1][1],r_soln_a[k][1][2]));
		refblock2.push_back(XYZ(r_soln_c[k][1][0],r_soln_c[k][1][1],r_soln_c[k][1][2]));
		refblock2.push_back(XYZ(r_soln_n[k][2][0],r_soln_n[k][2][1],r_soln_n[k][2][2]));
//		refblock2.push_back(XYZ(r_soln_a[k][2][0],r_soln_a[k][2][1],r_soln_a[k][2][2]));
		refblock2.push_back(fixedcrds[2]);
//		if(((block1[2]-refblock1[2]).squarednorm()+
//				(block1[3]-refblock1[3]).squarednorm()+
//				(block2[1]-refblock2[1]).squarednorm())< 0.00001) continue;
		QuatFit qf1,qf2;
		double rmsd21,rmsd22;
		rmsd21=qf1.setup(refblock1,block1);
		rmsd22=qf2.setup(refblock2,block2);
		solutions.push_back(LoopSolution());
		solutions.back().rt1=qf1.getRigidTransform();
		solutions.back().rt2=qf2.getRigidTransform();
		std::vector<NSPgeometry::XYZ> natomcrds=atomcrds;
		NSPgeometry::XYZ nc1,nn2, nca2,nca2_2,nc2,nn3,nca3,nc3;
		nc1=solutions.back().rt1.applytoCopy(atomcrds[C1]);
		nn2=solutions.back().rt1.applytoCopy(atomcrds[N2]);
		nca2=solutions.back().rt1.applytoCopy(atomcrds[CA2]);
		nc2=solutions.back().rt2.applytoCopy(atomcrds[C2]);
		nn3=solutions.back().rt2.applytoCopy(atomcrds[N3]);
		nca3=solutions.back().rt2.applytoCopy(atomcrds[CA3]);
		nc3=solutions.back().rt2.applytoCopy(atomcrds[C3]);
		double rcac_min=1000;
		double rcac_max=-1000;
		double rcac1=NSPgeometry::distance2(fixedcrds[1],nc1);
		rcac_min=rcac1<rcac_min?rcac1:rcac_min;
		rcac_max=rcac1>rcac_max?rcac1:rcac_max;
		double rcac2=NSPgeometry::distance2(nca2,nc2);
		rcac_min=rcac2<rcac_min?rcac2:rcac_min;
		rcac_max=rcac2>rcac_max?rcac2:rcac_max;
		bool accurate=true;
		if(rcac_min<1.48*1.48||rcac_max>1.58*1.58) accurate=false;

//		}
		double rnca3=NSPgeometry::distance2(nn3,fixedcrds[2]);
//		double rcaca3_old=NSPgeometry::distance2(atomcrds[CA3],fixedcrds[2]);
//		double rcaca3=NSPgeometry::distance2(nca3,fixedcrds[2]);
//		nca2_2=solutions.back().rt2.applytoCopy(atomcrds[CA2]);
//		double rdiff=NSPgeometry::distance2(nca2,nca2_2);
		if(rnca3<1.42*1.42||rnca3>1.52*1.52) accurate=false;
//			std::cout << "rnca3: " <<sqrt(rnca3)<<"\t"<< std::endl;
//					sqrt(rcaca3_old)<<"\t"<< sqrt(rcaca3)<<"\t"<<sqrt(rdiff)<<std::endl;
		if(!accurate) solutions.pop_back();

//		std::cout <<"--------"<<std::endl;
//		std::cout <<"coord n2: " <<nn2.toString() <<" " <<refblock1[2].toString() <<std::endl;
//		std::cout <<"coord ca2: " <<nca2.toString() <<" " <<refblock1[3].toString() <<std::endl;
//		std::cout <<"coord c2: " <<nc2.toString() <<" " <<refblock2[1].toString() <<std::endl;
//		std::cout <<"coord n3 " <<nn3.toString() <<" " <<refblock2[2].toString() <<std::endl;
//		std::cout <<"coord ca3: " <<nca3.toString() <<" " <<refblock2[3].toString() <<std::endl;
//		std::cout <<"coord c3: " <<nc3.toString() <<" " <<fixedcrds[3].toString() <<std::endl;
	}
	return solutions.size();
}


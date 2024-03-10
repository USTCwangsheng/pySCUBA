/*
 * quatfit.cpp
 *
 *  Created on: 2016年4月12日
 *      Author: hyliu
 */

#include "geometry/quatfit.h"
#include "geometry/rotation.h"


using namespace NSPgeometry;
double NSPgeometry::RMSD(const std::vector<NSPgeometry::XYZ> &refcrd,
		const std::vector<NSPgeometry::XYZ> &crd_before,
		std::vector<NSPgeometry::XYZ> &crd_after, const std::vector<double> weights){
	QuatFit qf;
	crd_after=crd_before;
	return sqrt(qf.fitting(refcrd,crd_after,weights));
}
void QuatFit::transform(std::vector<XYZ> & crd){
	for(auto & p:crd) p=p+shift_;
	Rotation r;
//	QuaternionCrd q=getquaternioncrd();
//	std::cout <<q.angle() <<std::endl;
//	std::cout <<shift_.toString() <<std::endl;
	r.init(getquaternioncrd(),ref_center_);

//	RigidTransform rt(q,ref_center_,r.applytoCopy(shift_));
	for(auto & p:crd) {
		r.apply(&p);
	}
}
double QuatFit::fitting(const std::vector<double> & ref_coord,
						std::vector<double> & coord,std::vector<double> w){
	std::vector<XYZ> refcrd,crd;
	XYZ::vectortoxyzs(ref_coord,refcrd);
	XYZ::vectortoxyzs(coord,crd);
	double rmsd2=fitting(refcrd,crd,w);
	coord.clear();
	XYZ::xyzstovector(crd,coord);
	return rmsd2;
}

double QuatFit::fitting(const std::vector<XYZ> & ref_coord,
						std::vector<XYZ> & coord,std::vector<double> w){
	double rmsd2=setup(ref_coord, coord, w);
	transform(coord);
	return rmsd2;
}

QuaternionCrd QuatFit::getquaternioncrd(){
	return QuaternionCrd({m_evec[0][0],m_evec[1][0],m_evec[2][0],m_evec[3][0]});
}

double QuatFit::setup(const std::vector<double> &ref_coord,const std::vector<double> & coord,
		std::vector<double> w){
	std::vector<XYZ> refcrd, crd;
	XYZ::vectortoxyzs(ref_coord,refcrd);
	XYZ::vectortoxyzs(coord,crd);
	return setup(refcrd,crd,w);
}

#include <lapacke.h>
//extern "C" lapack_int LAPACKE_dsyev( int matrix_layout, char jobz, char uplo, lapack_int n,
//                          double* a, lapack_int lda, double* w );
//lapack_int LAPACKE_dsyev( int matrix_layout, char jobz, char uplo, lapack_int n,
//                          double* a, lapack_int lda, double* w );
double QuatFit::setup( std::vector<XYZ> ref_coord,
						std::vector<XYZ> coord,std::vector<double> w) {
	// The rotation that fit coord onto refcoord
	    double a[3],b[3];
	    for (unsigned int i=0;i<4;i++){
	    	for(unsigned int j=0;j<=i;j++){
	    		m_matrix[i][j]=0.0;
	    	}
	    }
	    int np=ref_coord.size();
	    if(w.empty()) {
	    	for(int i=0; i<np; ++i) w.push_back(1.0);
	    }
	    ref_center_=center(ref_coord,w);
//	    std::cout <<ref_center_.toString() <<std::endl;

	    shift_=-1.0*center(coord,w);

//	    std::cout <<shift_.toString() <<std::endl;
	    for(int i=0;i<np; ++i) {
	    	ref_coord[i] = ref_coord[i]-ref_center_;
	    	coord[i]=coord[i]+shift_;
	    }
	    shift_=shift_+ref_center_;

	    m_wtot=0.0;
	    int n=w.size();
	    for (int i=0; i<n;i++ ){
	    	m_wtot += w[i];
	    	for (int j=0; j<3;j++) {
	    		a[j]= ref_coord[i][j]+coord[i][j];
	    		b[j]= ref_coord[i][j]-coord[i][j];
	    	}
	    	m_matrix[0][0]+=w[i]*(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
	    	m_matrix[1][0]+=w[i]*(a[2]*b[1]-a[1]*b[2]);
	    	m_matrix[2][0]+=w[i]*(-a[2]*b[0]+a[0]*b[2]);
	    	m_matrix[3][0]+=w[i]*(a[1]*b[0]-a[0]*b[1]);

	    	m_matrix[1][1]+=w[i]*(b[0]*b[0]+a[2]*a[2]+a[1]*a[1]);
	    	m_matrix[2][1]+=w[i]*(b[1]*b[0]-a[0]*a[1]);
	    	m_matrix[3][1]+=w[i]*(b[2]*b[0]-a[0]*a[2]);

	    	m_matrix[2][2]+=w[i]*(b[1]*b[1]+a[2]*a[2]+a[0]*a[0]);
	    	m_matrix[3][2]+=w[i]*(b[2]*b[1]-a[1]*a[2]);

	    	m_matrix[3][3]+=w[i]*(b[2]*b[2]+a[1]*a[1]+a[0]*a[0]);
	    }
	    double A[16];
	    for (unsigned int i=0;i<4;i++){
	    	for(unsigned int j=0;j<=i;j++){
	    		double x=m_matrix[i][j]/m_wtot;
	    		A[i*4+j]=x;
	    		A[j*4+i]=x;
	    	}
		}
	    lapack_int N{4};
	    lapack_int LDA{4};
	    lapack_int INFO;
	    INFO=LAPACKE_dsyev(LAPACK_COL_MAJOR,'V','U',N,A,LDA,m_eval);
	    for (unsigned int i=0;i<4;i++){
	    	for(unsigned int j=0;j<4;j++){
	    		m_evec[j][i]=A[i*4+j];
	    	}
		}
	    return m_eval[0];
}


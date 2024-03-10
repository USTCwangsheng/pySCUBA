/*
 * pca.cpp
 *
 *  Created on: 2017年6月14日
 *      Author: hyliu
 */




#include "dstl/pca.h"
#include <cmath>

using namespace NSPdstl;
PCA::PCA(int ndim, bool angular) :
		ndim_(ndim), angular_(angular) {
	center_.resize(ndim, 0);
	covarmatrix_.resize(ndim, ndim, 0);
	components_.resize(ndim, ndim, 0);
	eigenvalues_.resize(ndim, 0);
}
std::vector<double> PCA::getPCAcrds(const std::vector<double> &point,int ndim) const {
	assert(ndim<=ndim_);
	std::vector<double> res;
	if(ndim==0) ndim=ndim_;
	for(int i=0; i<ndim; ++i) {
		const std::vector<double> & cpnti=getcomponent(i);
		double p{0};
		for(int j=0; j<ndim_;++j) {
			double delta= diff(point[j],center_[j]);
			p +=delta*cpnti[j];
		}
		res.push_back(p);
	}
	return res;
}
void PCA::calccenter(const std::vector<std::vector<double>> &data) {
	double datasize = 1.0 / (double) data.size();
	center_.resize(ndim_, 0.0);
	if (!angular_) {
		for (auto & d : data) {
			for (int j = 0; j < ndim_; j++) {
				center_[j] += d.at(j);
			}
		}
		for (auto &c : center_)
			c = c * datasize;
	} else {
		std::vector<double> xav(ndim_, 0.0);
		std::vector<double> yav(ndim_, 0.0);
		double deg = 3.1415926535897932384626 / 180.0;
		for (auto &d : data) {
			for (int j = 0; j < ndim_; j++) {
				double ang = d.at(j) * deg;
				xav[j] += cos(ang);
				yav[j] += sin(ang);
			}
		}
		for (int j = 0; j < ndim_; j++) {
			double x = xav[j] * datasize;
			double y = yav[j] * datasize;
			if (x > -1.e-10 && x < 1.e-10) {
				if (y >= 0)
					center_[j] = 90.0;
				else
					center_[j] = -90.0;
			} else {
				double ang = atan(y / x) / deg;
				if (x < 0)
					ang = ang + 180.0;
				if (ang > 180.0)
					ang -= 360.0;
				center_[j] = ang;
			}
		}
	}
}
void PCA::calccovarmatrix(const std::vector<double> &center,const std::vector<std::vector<double>> &data){
	center_=center;
	MyMatrix sum2dev;
	sum2dev.resize(ndim_,ndim_,0.0);
	for(auto &d:data) {
		for(int i=0; i<ndim_; ++i) {
			double devi=diff(d.at(i),center.at(i));
			sum2dev(i,i) += devi*devi;
			for(int j=0;j<i; ++j) {
				sum2dev(i,j) += devi*diff(d.at(j),center.at(j));
			}
		}
	}
	double ndata_i=1.0/(double) data.size();
	for(int i=0;i<ndim_; ++i) {
		for(int j=0; j<=i; ++j) {
			covarmatrix_(i,j)=sum2dev(i,j)*ndata_i;
			if(j !=i) covarmatrix_(j,i)=covarmatrix_(i,j);
		}
	}
}

#include <lapacke.h>
MyMatrix NSPdstl::diagonizematrix(const MyMatrix &mat,std::vector<double> *eigenvalues){
	int ndim=mat.M;
    double *A=new double[ndim*ndim];
    for (unsigned int i=0;i<ndim;i++){
    	for(unsigned int j=0;j<ndim;j++){
    		A[i*ndim+j]=mat(i,j);
    	}
	}
    lapack_int N=ndim;
    lapack_int LDA=ndim;
    lapack_int INFO;
    double *eval=new double[ndim];
    INFO=LAPACKE_dsyev(LAPACK_COL_MAJOR,'V','U',N,A,LDA,eval);
    eigenvalues->resize(ndim,0);
    for(unsigned int i=0;i<ndim;++i){
    	eigenvalues->at(i)=eval[ndim-i-1];
    }
    MyMatrix res;
    res.resize(ndim,ndim,0);
    for (unsigned int i=0;i<ndim;i++){
    	std::vector<double> v;
    	for(unsigned int j=0;j<ndim;j++){
    		v.push_back(A[i*ndim+j]);
    	}
    	res.setcol(ndim-i-1,v);
	}
   delete[] A;
   delete[] eval;
   return res;
}


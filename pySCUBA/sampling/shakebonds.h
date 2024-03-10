/*
 * shakebonds.h
 *
 *  Created on: 2019年12月16日
 *      Author: hyiu
 */

#ifndef SHAKEBONDS_H_
#define SHAKEBONDS_H_
#include "geometry/calculators.h"
#include <memory>
#include <functional>
namespace  NSPintrct{
	class IntrctMol;
};
namespace NSPsampling{
class ShakeBonds {
public:
	struct Bond{
		int i;
		int j;
		double b02;
		Bond(int a1,int a2, double b): i(a1),j(a2),b02(b*b){;}
	};
	bool shake(const std::vector<double> & masses,const std::vector<double> &crdref, std::vector<double> &crd) const;
	std::function<bool (const std::vector<double> &,const std::vector<double> &, std::vector<double> &)>
			shakefunction() {
		return std::bind(&ShakeBonds::shake, *this, std::placeholders::_1,
				std::placeholders::_2,std::placeholders::_3);
	}
	void addbond(int a1,int a2, double b){
		bonds.push_back(Bond(a1,a2,b));
	}
	int & ndim() {return ndim_;}
	const int & ndim() const {return ndim_;}
	int nbonds() const {return bonds.size();}
	void seton(bool val=true){ison_=val;}
private:
	std::vector<Bond> bonds;
	bool ison_{true};
	int ndim_{3};
	std::vector<double> rij(const std::vector<double> & crd, int i, int j) const{
		int iidx=i*ndim_;
		int jidx=j*ndim_;
		std::vector<double> diff;
		for(int d=0;d<ndim_;++d){
			diff.push_back(crd[jidx+d]-crd[iidx+d]);
		}
		return diff;
	}
	double norm2(const std::vector<double> &rij) const {
		double r2=0.0;
		for(int d=0;d<ndim_;++d){
			r2+=rij[d]*rij[d];
		}
		return r2;
	}
	double dot(const std::vector<double> & ra, const std::vector<double> &rb) const {
		double res=0.0;
		for(int d=0;d<ndim_;++d){
			res +=ra[d]*rb[d];
		}
		return res;
	}
};

std::shared_ptr<ShakeBonds> make_shakebonds(NSPintrct::IntrctMol &mol);
}


#endif /* SHAKEBONDS_H_ */

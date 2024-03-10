/*
 * stochasticdynamics.h
 *
 *  Created on: 2019年12月16日
 *      Author: hyiu
 */

#ifndef STOCHASTICDYNAMICS_H_
#define STOCHASTICDYNAMICS_H_
#include "dstl/randomengine.h"
#include <vector>
#include <cmath>
#include <functional>
#include <memory>
#include <cassert>
#define MASS_MAX 0.1*DBL_MAX
namespace NSPintrct{
  class IntrctMol;
}
namespace NSPsampling {
struct SDRunPara;
class StochasticDynamics {
public:
	struct State {
		std::vector<double> crd;
		std::vector<double> vel;
		std::vector<double> sdintegral;
		void scaletemp(double scale){
			double scalert=sqrt(scale);
			for(auto & v:vel) v*=scalert;
			for(auto &si:sdintegral) si*=scalert;
		}
		void scaletemp(double scale, const std::vector<int> atoms, int ndim=3){
			double scalert=sqrt(scale);
			for(int a:atoms){
				for(int m=0;m<ndim;++m){
					int id=ndim*a+m;
					vel[id] *=scalert;
					sdintegral[id] *=scalert;
				}
			}
		}
	};
	struct SDStruct {
		double c1, c2, c3, c4, c5, c6, c7, c8, c9;
		SDStruct(double gamma, double timestep);
	};
	typedef std::function<
			bool(const std::vector<double> &, const std::vector<double> &,
					std::vector<double> &)> ShakeFunction;
	StochasticDynamics() {
		;
	}
	void init(const SDRunPara &para, const NSPintrct::IntrctMol &imol);
	void init(const std::vector<double> &masses,
			const std::vector<double> & gammas, double timestep, double kbT,
			int ndim = 3);
	void init(const std::vector<double> &masses,
			const std::vector<double> & gammas, double timestep,
			const std::vector<std::vector<int>> &atomgroups,const std::vector<double> & kbT,
			int ndim = 3);
	bool shakestate(State & state, const ShakeFunction &shake, bool shakecrd,
			const std::vector<double> *vmasses=nullptr);
	std::shared_ptr<State> make_initstate(const std::vector<double> &initcrd,NSPdstl::RandomEngine<> &rng,const std::vector<double> *vmasses=nullptr);
	bool leapstep(State &stateold, State &statenew,NSPdstl::RandomEngine<> &rng,
			const std::vector<double> & forces, const ShakeFunction & shake,
			const std::vector<double> *vmasses=nullptr) const;
	double & timestep() {
		return timestep_;
	}
	double ekin(const State &state) const {
		double res = 0.0;
		for (int i = 0; i < ndof_; ++i) {
			double v = state.vel[i];
			res += 0.5 * masses_[i] * v * v;
		}
		return res;
	}
	std::vector<double> & kbT() {
		return kbT_;
	}
	const std::vector<double> &kbT() const{
		return kbT_;
	}
	void assigntemperatures(const std::vector<std::vector<int>> & atomgroups,
			const std::vector<double> &newkbT,int ndim=3){
		assert(atomgroups.size()==newkbT.size());
		kbT_.resize(ndof_,1.0);
		int nassigned=0;
		for(int g=0;g<atomgroups.size();++g){
			for(int a:atomgroups[g])
				for(int m=0;m<ndim;++m)kbT_[ndim*a+m]=newkbT[g];
			nassigned+=atomgroups[g].size()*ndim;
		}
		assert(nassigned==ndof_);
	}
	void scaletemperatures(double factor,
			const std::vector<int> &atoms=std::vector<int>(),int ndim=3){
		if(atoms.empty()){
			for(double &t:kbT_) t *=factor;
		}
		for(int a:atoms)
			for(int m=0;m<ndim;++m) kbT_[ndim*a+m] *=factor;
	}
	const std::vector<double> &masses() const {return masses_;}
private:
	int ndof_ { 0 };
//	double kbT_ { 1.0 };
	double timestep_ { 0.002 };
	std::vector<double> kbT_;
	std::vector<double> masses_;
	std::vector<double> gammas_;
	std::vector<SDStruct> sdstructs_;
};
}


#endif /* STOCHASTICDYNAMICS_H_ */

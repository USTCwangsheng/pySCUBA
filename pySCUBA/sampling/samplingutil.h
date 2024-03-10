/*
 * samplingutil.h
 *
 *  Created on: 2019年12月24日
 *      Author: hyiu
 */

#ifndef SAMPLINGUTIL_H_
#define SAMPLINGUTIL_H_
#include "iblock/intrctmol.h"
#include "sampling/molsdrun.h"
#include <cassert>
#include <cmath>
namespace NSPsampling {
/**
 * stop criteria in a sampling process
 *
 * The criteria are met when the number of steps reached maxcycles_, or
 * if the mode is VAR, variations of time-averaged energy fluctuation below enevarcut_
 * if the model is MINPERSIST, the persistence steps of the minimum energy sampled so far
 */
struct StopJudger {
	/**
	 * keep tracks of minimum energy and its persistence steps
	 */
	struct MinEneKeeper {
		double minene { 1.e16 };
		int pcycles { 0 }; //<* number of cycles passed since last update
		bool update(double e) {
			if (e > minene) {
				++pcycles;
				return false;
			}
			minene = e;
			pcycles = 0;
			return true;
		}
	};
	/**
	 * determine and keep time-averaged energy and energy fluctuation
	 */
	struct EneAverager {
		double decay;
		double e_average { 0.0 };
		double e2_average { 0.0 };
		double decay_acc { 1.0 };
		void init(double d) {
				decay=d;
				e_average=0.0;
				e2_average=0.0;
				decay_acc=1.0;
		}
		void accumulate(double ene) {
			e_average = decay * e_average +(1.0-decay)* ene;
			e2_average = decay * e2_average + (1.0-decay)*ene*ene;
			if (decay_acc > 1.e-10)
				decay_acc *= decay;
		}
		double fluct2() const {
			if (decay_acc < 1.e-3)
				return e2_average - e_average * e_average;
			return 1000000.0;
		}
	};
	enum JUDGEMODE{VAR,MINPERSIST};
	StopJudger(){;}
	StopJudger(int mode,int maxcyc,double d1,double d2){
		assert(mode==VAR || mode==MINPERSIST);
		mode_=mode;
		maxcycles_=maxcyc;
		if(mode_==VAR){
			averager_.init(d1);
			enevarcut_=d2*d2;
		} else {
			 eminpersiscycles_=(int) d1;
		}
	}
	bool operator()(int ncyc, double e){
		if(stopbycycles(ncyc)) return true;
		if(mode_==VAR) return stopbyenevar(e);
		else return stopbyeminpersist(e);
	}
	bool stopbycycles(int ncyc){
		if(ncyc>=maxcycles_) return true;
		return false;
	}
	bool stopbyeminpersist(double e){
		 minkeeper_.update(e);
		 if(minkeeper_.pcycles >=eminpersiscycles_) return true;
		 return false;
	}
	int eminpersist() const {
		return minkeeper_.pcycles;
	}
	bool stopbyenevar(double e){
		averager_.accumulate(e);
		if(averager_.fluct2()<=enevarcut_) return true;
		return false;
	}
	double enevar() const{
		return sqrt(averager_.fluct2());
	}
private:
	EneAverager averager_;
	MinEneKeeper minkeeper_;
	int mode_{MINPERSIST};
	int maxcycles_{-1}; //<* maximum number of cycles;
	double enedecay_{1.0}; //<*per-step decay for estimating energy variations by time averaging;
	double enevarcut_{0.0}; //<* timer-averaged energy variation-based  stop criterion
	int eminpersiscycles_{-1}; //*<prsistence steps of minimum energy-based stop criterion.
};
class SDCallBackOptSD{
public:
	mutable StopJudger judger;
	int printsteps{10};
	bool verbose{true};
	mutable std::array<double,NSPintrct::IntrctBlck::ENESIZE> energies;
    template <typename PARA>
	SDCallBackOptSD(const PARA &para){
		judger=para.mksdjudger();
		printsteps=para.sdprintsteps;
		verbose=para.verbose;
	}
	bool operator() (MolSDRun &run) const{
		double epot_tot=run.imol()->sum_energies(energies);
		if((run.nstepsrun()+1)%printsteps==0){
				double temp=run.temperature();
				std::cout<<"Step "<< run.nstepsrun()
					<<"Temperature= "<<temp
					<<" Pot_tot= "<<epot_tot
					<<" Restraint_ene: ";
						for(auto &r:run.restraints()) {
							std::cout << " "<<r->ene();
						}
					std::cout<<std::endl;
			}
		bool tostop= judger(run.nstepsrun()+1,epot_tot);
		if(verbose){
			std::cout <<"SD steps " <<run.nstepsrun()+1 << " ENE variation: "<< judger.enevar()<<std::endl;
		}
		return tostop;
	}
};

}

#endif /* SAMPLINGUTIL_H_ */

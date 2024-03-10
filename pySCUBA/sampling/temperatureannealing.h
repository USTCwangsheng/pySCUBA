/*
 * temperatureannealing.h
 *
 *  Created on: 2019年12月23日
 *      Author: hyiu
 */

#ifndef TEMPERATUREANNEALING_H_
#define TEMPERATUREANNEALING_H_
#include <cassert>
namespace NSPsampling{
/**
 * data defining how bath temperature changes in simulated annealing SD
 */
struct TemperatureAnnealing{
	enum PHASE{HIGHT,CHANGING,LOWT };
	double t_high{0.0}; //*< high temperature
	double t_low{0.0}; //*< low temperature
	int nsteps_cycle{0}; //*< total steps of a complete annealing cycle
	int nsteps_high{0}; //*< number of steps in high temperature
	int nsteps_drop{0}; //*< number of steps to change gradually from high to low
	int annealinggrp{-1};
	TemperatureAnnealing(){;}
	TemperatureAnnealing(double th,double tl, double cycle,double sh,double sd,int agrp=-1):
	t_high(th),t_low(tl),nsteps_cycle(cycle),nsteps_high(sh),nsteps_drop(sd),annealinggrp(agrp){
		assert(sh+sd <= cycle);
	}
	/**
	 * determine the current bath temperature
	 */
	double temperature(int step) const {
		step=(step-1)%nsteps_cycle;
		if(step<nsteps_high) return t_high;
		else if(step>=(nsteps_high+nsteps_drop)) return t_low;
		else return t_high - (t_high-t_low)*(double)(step-nsteps_high)/(double)(nsteps_drop);
	}
	PHASE phase(int step){
		step=(step-1)%nsteps_cycle;
		if(step <nsteps_high) return HIGHT;
		else if(step>=(nsteps_high+nsteps_drop)) return LOWT;
		else return CHANGING;
	}
};
}


#endif /* TEMPERATUREANNEALING_H_ */

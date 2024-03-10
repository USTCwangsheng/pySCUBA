/*
 * stochasticdynamics.cpp
 *
 *  Created on: 2019年12月16日
 *      Author: hyiu
 */
#include "sampling/stochasticdynamics.h"
#include "sampling/molsdrun.h"
#include "sampling/shakebonds.h"
#include <cfloat>
#define MASS_MAX 0.1*DBL_MAX
using namespace NSPsampling;
void StochasticDynamics::init(const SDRunPara &para, const NSPintrct::IntrctMol & imol){
	std::vector<double> masses;
	std::vector<double> gammas;
	double gamma=para.gamma;
	double atommass=para.atommass;
	masses.assign(imol.natoms(),atommass);
	gammas.assign(imol.natoms(),gamma);
	if(!para.atomgroups.empty()){
		init(masses, gammas, para.timestep, para.atomgroups, para.temperatures,3);
	} else {
		init(masses,gammas,para.timestep,para.temperatures[0],3);
	}
}
StochasticDynamics::SDStruct::SDStruct(double gamma, double timestep) {
	double gdt = gamma * timestep;
	double gdth = 0.5 * gdt;
	if (gdt > 0.05) {
		const double emdth = exp(-gdth);
		const double epdth = exp(+gdth);
		const double emdt = emdth * emdth;
		const double epdt = epdth * epdth;
		const double omdt = 1.0 - emdt;
		const double cdth = gdt - 3.0 + 4.0 * emdth - emdt;
		const double ddth = 2.0 - epdth - emdth;
		const double bpdth = gdt * (epdt - 1.0)
				- 4.0 * (epdth - 1.0) * (epdth - 1.0);
		const double bmdth = gdt * (1.0 - emdt)
				- 4.0 * (emdth - 1.0) * (emdth - 1.0);
		c1 = emdt;
		c2 = omdt / gdt;
		c3 = sqrt(fabs(omdt));
		c4 = sqrt(fabs(bpdth / cdth));
		c5 = gamma * ddth / cdth;
		c6 = (epdth - emdth) / gdt;
		c7 = sqrt(fabs(cdth)) / gamma;
		c8 = sqrt(fabs(bmdth / omdt)) / gamma;
		c9 = -ddth / (gamma * omdt);
	} else {
		const double gdth2 = gdth * gdth;
		const double gdth3 = gdth2 * gdth;
		const double gdth4 = gdth2 * gdth2;
		const double gdth5 = gdth2 * gdth3;
		const double gdth6 = gdth3 * gdth3;
		const double gdth7 = gdth4 * gdth3;
		c1 = exp(-gdt);

		c2 = 1.0 - gdth + gdth2 * 2.0 / 3.0 - gdth3 / 3.0 + gdth4 * 2.0 / 15.0
				- gdth5 * 2.0 / 45.0 + gdth6 * 4.0 / 315.0;
		c3 = sqrt(c2 * 2.0 * gdth);

		c4 = sqrt(
				fabs(
						gdth / 2.0 + gdth2 * 7.0 / 8.0 + gdth3 * 367.0 / 480.0
								+ gdth4 * 857.0 / 1920.0
								+ gdth5 * 52813.0 / 268800.0
								+ gdth6 * 224881.0 / 3225600.0
								+ gdth7 * 1341523.0 / 64512000.0));

		c5 = -2.0 / timestep
				* (1.5 + gdth * 9.0 / 8.0 + gdth2 * 71.0 / 160.0
						+ gdth3 * 81.0 / 640.0 + gdth4 * 7807.0 / 268800.0
						+ gdth5 * 1971.0 / 358400.0
						+ gdth6 * 56417.0 / 64512000.0);

		c6 = 1.0 + gdth2 / 6.0 + gdth4 / 10.0 + gdth6 / 5040.0;

		c7 = timestep * 0.5
				* sqrt(
						fabs(
								gdth * 2.0 / 3.0 - gdth2 / 2.0
										+ gdth3 * 7.0 / 30.0 - gdth4 / 12.0
										+ gdth5 * 31.0 / 1260.0 - gdth6 / 160.0
										+ gdth7 * 127.0 / 90720.0));

		c8 = timestep * 0.5
				* sqrt(
						fabs(
								gdth / 6.0 - gdth3 / 60.0
										+ gdth5 * 17.0 / 10080.0
										- gdth7 * 31.0 / 181440.0));

		c9 = timestep * 0.5
				* (0.5 + gdth / 2.0 + gdth2 * 5.0 / 24.0 + gdth3 / 24.0
						+ gdth4 / 240.0 + gdth5 / 720.0 + gdth6 * 5.0 / 8064.0);
	}
}
void StochasticDynamics::init(const std::vector<double> & masses,
		const std::vector<double> & gammas, double timestep, double kbT,
		int ndim) {
	assert(gammas.size() == masses.size());
	ndof_ = masses.size() * ndim;
	timestep_ = timestep;
	kbT_.assign(ndof_,kbT);
	masses_.clear();
	gammas_.clear();
	for (int i = 0; i < masses.size(); ++i) {
		for (int m = 0; m < ndim; ++m) {
			masses_.push_back(masses[i]);
			gammas_.push_back(gammas[i]);
		}
	}
	sdstructs_.clear();
	for (int i = 0; i < ndof_; ++i) {
		sdstructs_.push_back(SDStruct(gammas_[i], timestep_));
	}
}
void StochasticDynamics::init(const std::vector<double> & masses,
		const std::vector<double> & gammas, double timestep,
		const std::vector<std::vector<int>> &atomgroups,
		const std::vector<double>& kbT,
		int ndim) {
	init(masses,gammas,timestep,1.0,ndim);
	assigntemperatures(atomgroups,kbT);
}
std::shared_ptr<StochasticDynamics::State> StochasticDynamics::make_initstate(
		const std::vector<double> &initcrd, NSPdstl::RandomEngine<> &rng,
		const std::vector<double> *vmasses) {
	assert(initcrd.size() == ndof_);
	auto state = std::shared_ptr < State > (new State);
	std::vector<double> &crd = state->crd;
	std::vector<double> &vel = state->vel;
	std::vector<double> &sdintegral = state->sdintegral;
	const std::vector<double> *m=&masses_;
	if(vmasses) m=vmasses;
	const std::vector<double> & mss=*m;
	crd.resize(ndof_);
	std::copy(initcrd.begin(), initcrd.end(), crd.begin());
	for (int i = 0; i < ndof_; ++i) {
		double sd = sqrt(kbT_[i] / mss[i]);
		vel.push_back(rng.randomnormal(0.0, sd));
		sd = sdstructs_[i].c7 * (kbT_[i] / mss[i]);
		sdintegral.push_back(rng.randomnormal(0.0, sd));
	}
	return state;
}
bool StochasticDynamics::shakestate(State &state, const ShakeFunction & shake,
		bool shakecrd,const std::vector<double> *vmasses) {
	std::vector<double> &crd = state.crd;
	std::vector<double> &vel = state.vel;
	std::vector<double> crd_old = crd;
	const std::vector<double> *m=&masses_;
	if(vmasses) m=vmasses;
	const std::vector<double> &mss=*m;
	if (shakecrd)
//		if (!shake(masses_, crd_old, crd)) return false;
		if (!shake(*m, crd_old, crd)) return false;
//	std::vector<double> &vel=state_->vel;
//	std::vector<double> &sdintegral=state_->sdintegral;
//	auto  & rng=NSPdstl::RandomEngine<>::getinstance();
	for (int i = 0; i < ndof_; i++) {
		if(mss[i]==MASS_MAX) continue;
		double kToverM = sqrt(kbT_[i] / mss[i]);
		SDStruct &sdstruct = sdstructs_[i];
		double cf = timestep_ * sdstruct.c2 / mss[i];
		double sd1 = sdstruct.c3 * kToverM;
		double sd2 = kToverM * sdstruct.c4;
		crd_old[i] = crd[i] - vel[i] * timestep_ * sdstruct.c6;
	}
//	if (!shake(masses_, crd, crd_old)) return false;
	if (!shake(*m, crd,crd_old)) return false;
	for (int i = 0; i < ndof_; ++i) {
		if(mss[i]==MASS_MAX) continue;
		SDStruct &sdstruct = sdstructs_[i];
		vel[i] = (crd[i] - crd_old[i]) / (timestep_ * sdstruct.c6);
	}
	return true;
}

bool StochasticDynamics::leapstep(State &stateold, State &statenew,NSPdstl::RandomEngine<> &rng,
		const std::vector<double> & forces, const ShakeFunction & shake,
		const std::vector<double> *vmasses) const {
	std::vector<double> &crd_old = stateold.crd;
	std::vector<double> &vel_old = stateold.vel;
	std::vector<double> &sdintegral_old = stateold.sdintegral;
	std::vector<double> &crd = statenew.crd;
	std::vector<double> &vel = statenew.vel;
	std::vector<double> &sdintegral = statenew.sdintegral;
	const std::vector<double> *m=&masses_;
	if(vmasses) m=vmasses;
	const std::vector<double> &mss=*m;
	for (int i = 0; i < ndof_; i++) {
		if(mss[i]==MASS_MAX) continue;
		double kToverM = sqrt(kbT_[i] / mss[i]);
		const SDStruct &sdstruct = sdstructs_[i];
		double cf = timestep_ * sdstruct.c2 / mss[i];
		double sd1 = sdstruct.c3 * kToverM;
		double sd2 = kToverM * sdstruct.c4;
		double vrand1 = rng.randomnormal(0.0, sd1);
		double vrand2 = rng.randomnormal(0.0, sd2);
		double svh = sdintegral_old[i] * sdstruct.c5 + vrand2;
		sdintegral[i] = vrand1;
		//2.11.2.2, 2.11.2.3
		vel[i] = (vel_old[i] - svh) * sdstruct.c1 + forces[i] * cf + vrand1;
		crd[i] = crd_old[i] + vel[i] * timestep_ * sdstruct.c6;
	}
//	if (!shake(masses_, crd_old, crd)) {
	if (!shake(mss, crd_old, crd)) {
		return false;
	}
	for (int i = 0; i < ndof_; ++i) {
		if(mss[i]==MASS_MAX) continue;
		double kToverM = sqrt(kbT_[i] / mss[i]);
		const SDStruct &sdstruct = sdstructs_[i];
		vel[i] = (crd[i] - crd_old[i]) / (timestep_ * sdstruct.c6);
		double sd3 = kToverM * sdstruct.c7;
		double sd4 = kToverM * sdstruct.c8;
		double vrand3 = rng.randomnormal(0.0, sd3);
		double vrand4 = rng.randomnormal(0.0, sd4);
		double sxh = sdintegral[i] * sdstruct.c9 + vrand4;
		sdintegral[i] = vrand3;
		crd[i] += vrand3 - sxh;
	}
//	if (!shake(masses_, crd_old, crd))
	if (!shake(mss, crd_old, crd))
		return false;
	return true;
}


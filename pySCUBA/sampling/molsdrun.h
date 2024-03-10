/*
 * molsdrun.h
 *
 *  Created on: 2019年12月16日
 *      Author: hyiu
 */

#ifndef MOLSDRUN_H_
#define MOLSDRUN_H_
#include "sampling/stochasticdynamics.h"
#include "iblock/intrctmol.h"
#include "sampling/restraints.h"
#include "sampling/sdrunpara.h"
#include "iblock/molsystmpara.h"
#include "sampling/shakebonds.h"
#include "dstl/randomengine.h"
#include "sampling/temperatureannealing.h"
namespace NSPsampling {
class MolSDRun;

/**
 * callback function after finishing a SD integration step
 */
class EmptySDCallBack {
public:
	bool operator()(MolSDRun &run) const {
		return false;
	}
};

/**
 * object driving a stochastic simulation of a molecular system.
 *
 * By using smart pointers to the Molsystm object
 * as data member, it can effectively hold this object.
 */
class MolSDRun {
public:
	MolSDRun() {
		;
	}

	template<typename SEED>
	void initrandomengine(SEED seed) {
		rng_ = std::shared_ptr<NSPdstl::RandomEngine<>>(
				new NSPdstl::RandomEngine<>);
		rng_->init(seed);
	}
	std::shared_ptr<StochasticDynamics> sd() {
		return sd_;
	}
	std::shared_ptr<NSPintrct::IntrctMol> imol() {
		return imol_;
	}
	std::shared_ptr<ShakeBonds> shakebds() {
		return shakebds_;
	}
	/**
	 * prepare to run SD on a pre-constructed molecular system.
	 */
	bool init(const SDRunPara &para, std::shared_ptr<NSPintrct::IntrctMol> imol,
			bool shakeinit = true);

	void addrestraint(std::shared_ptr<ResTraint> &restr) {
		restraints_.push_back(restr);
	}
	/**
	 * run the SD for a number of steps.
	 * After each step the supplied callback function is called.
	 */
	template<typename CALLBACK = EmptySDCallBack>
	bool runsteps(int nsteps, const CALLBACK &stepcallback = CALLBACK()) {
		std::vector<NSPgeometry::XYZ> dedx;
		for (int i = 0; i < nsteps; ++i) {
			if (nstepsrun_ == 0) {
				imol_->my_intrctpara->phipsicodes_updateonly = false;
				imol_->my_intrctpara->sscodes_updateonly = false;
				imol_->my_intrctpara->calc_nblists = true;
				imol_->my_intrctpara->calc_sscodes = true;
			} else {
				imol_->my_intrctpara->phipsicodes_updateonly = true;
				imol_->my_intrctpara->sscodes_updateonly = true;
				imol_->my_intrctpara->calc_nblists = false;
				if (!this->restrainhelices_)
					imol_->my_intrctpara->calc_sscodes = false;
				if (nstepsrun_ % nblsteps_ == 0)
					imol_->my_intrctpara->calc_nblists = true;
				if (nstepsrun_ % sscodesteps_ == 0)
					imol_->my_intrctpara->calc_sscodes = true;
			}
			std::vector<NSPgeometry::XYZ> xyz = NSPgeometry::doublevtoXYZv(
					state_->crd);
			imol_->forces_all(xyz, dedx);
			/*//debug test forces
			 std::array<double,NSPintrct::IntrctBlck::ENESIZE> energies;
			 std::vector<NSPgeometry::XYZ> dedxtmp;
			 for(int i=0;i<xyz.size();++i){
			 for (int m=0;m<3;++m){
			 if(dedx[i][m]==0) continue;
			 xyz[i][m] +=0.0001;
			 imol_->forces_all(xyz,dedxtmp);
			 double ep=imol_->sum_energies(energies);
			 xyz[i][m] -=0.0002;
			 imol_->forces_all(xyz,dedxtmp);
			 double em=imol_->sum_energies(energies);
			 std::cout << i <<" " <<m <<" " << (ep-em)/0.0002 <<"  "<<dedx[i][m]<<std::endl;
			 xyz[i][m] +=0.0001;
			 }
			 }
			 //end debug test forces */
			/*			if(this->restrainhelices_) {
			 imol_->my_intrctpara->calc_sscodes=true;
			 if(this->nstepsrun_==0){
			 restraints_.push_back(std::shared_ptr<HelixRestraint>(new HelixRestraint(imol_.get(),
			 initpara_.flankingmcoff)));
			 }
			 }*/

			for (auto &f : dedx)
				f = -f;

//in case of annealing and rg restraint, set krg=0 at high and dropping temperature phase
// to allow backbone to expand
			double krg_orig;
			int resswitch = 0;

			if (initpara_.doannealing) {
				TemperatureAnnealing ta(initpara_.annealingscheme[0],
						initpara_.annealingscheme[1],
						initpara_.annealingscheme[2],
						initpara_.annealingscheme[3],
						initpara_.annealingscheme[4]);
				int phase = ta.phase(nstepsrun_ + 1);
				//for (auto &r : restraints_) {
				//	auto x = dynamic_cast<RgRestraint *>(r.get());
				//	if (x)
				//		krg_orig = x->kres();
				//}
				//double krg = 0.0;
				//if (phase == TemperatureAnnealing::LOWT) {
				//	krg = krg_orig;
				//}
				if (phase != TemperatureAnnealing::LOWT) {
					resswitch = 1;
				}
				//for (auto &r : restraints_) {
				//	auto x = dynamic_cast<RgRestraint *>(r.get());
				//	if (x)
				//		x->kres() = krg;
				//}
			}
			for (auto &r : restraints_) {
				if (r->switchable_ == 1 && resswitch == 1) {
					r->energy(xyz, &dedx, resswitch);
				}
				else r->energy(xyz, &dedx);
			}
			//if (initpara_.doannealing) {
			//	for (auto &r : restraints_) {
			//		auto x = dynamic_cast<RgRestraint *>(r.get());
			//		if (x)
			//			x->kres() = krg_orig;
			//	}
			//}
/*			if (initpara_.avoidlongstrand) {
				if (nstepsrun_ % initpara_.ncheckstrandsteps == 0)
					chngstrandflankingweights();
			}*/
			std::vector<double> forces = NSPgeometry::XYZvtodoublev(dedx);
			bool done = stepcallback(*this);
			if (done)
				return true;
			if (!sd_->leapstep(*state_, *buffstate_, *rng_, forces,
					shakebds_->shakefunction(), &vmasses_)) {
				std::cout << "Shake failure took place" << std::endl;
				return false;
			}
			auto temp = state_;
			state_ = buffstate_;
			buffstate_ = temp;
			++nstepsrun_;
		}
		return true;
	}
	StochasticDynamics::State & state() {
		return *state_;
	}
	/**
	 * reinitialize the simulator with new coordinates and
	 * newrandomseed (if not zero). Reset nstepsrun_ to zero.
	 */
	void reinitstate(const std::vector<double> &newcrd, int newrandomseed = 0) {
		if (newrandomseed != 0)
			rng_->reseed(newrandomseed);
		state_ = sd_->make_initstate(newcrd, *rng_, &vmasses_);
		nstepsrun_ = 0;
	}
	const StochasticDynamics::State &state() const {
		return *state_;
	}
	int nstepsrun() const {
		return nstepsrun_;
	}

	/**
	 * calculate the current temperature from kinetic energies
	 */
	double temperature() const {
		double nfreedof;
		if (shakeon_)
			nfreedof = state_->crd.size() - nfixedcrd_ - shakebds_->nbonds();
		else
			nfreedof = state_->crd.size() - nfixedcrd_;
		return 2.0 * (sd_->ekin(*state_)) / (nfreedof * KBT);
	}

	bool & shakeon() {
		return shakeon_;
	}
	const bool & shakeon() const {
		return shakeon_;
	}
	/**
	 * every number of steps to recalculate neighbor list
	 */
	int &nblsteps() {
		return nblsteps_;
	}
	const int & nblsteps() const {
		return nblsteps_;
	}
	void chngintrctbyss();
	/**
	 * change the environment temperature
	 *
	 * atomic velocities are also scaled upon to ratio between new and old environment temperatures
	 */
	void changetemperature(double kbt_new) {
		double scale = kbt_new * KBT / sd_->kbT()[0];
		state_->scaletemp(scale);
		sd_->scaletemperatures(scale);
	}

	/**
	 * change environment temperature for a given group of atoms
	 */
	void changetemperature(double kbt_new, const std::vector<int> &atoms) {
		double scale = kbt_new * KBT / sd_->kbT()[3 * atoms[0]];
		state_->scaletemp(scale, atoms);
		sd_->scaletemperatures(scale, atoms, 3);
	}

	/**
	 * change the environmental temperature for atoms belonging to a temperature group
	 */
	void changetemperature(double kbt_new, int tgrp) {
		changetemperature(kbt_new, temperaturegroups_->at(tgrp));
	}
	std::shared_ptr<std::vector<std::vector<int>>> & temperaturegroups() {
		return temperaturegroups_;
	}
	const std::shared_ptr<std::vector<std::vector<int>>> & temperaturegroups() const {
		return temperaturegroups_;
	}
	std::vector<double> &bathtemperatures() {
		return bathtemperatures_;
	}
	const std::vector<std::shared_ptr<ResTraint> > & restraints() const {
		return restraints_;
	}
	void chngstrandflankingweights();
	const SDRunPara & initpara() const {
		return initpara_;
	}
private:
	std::shared_ptr<StochasticDynamics> sd_ { nullptr }; //*< the SD engine, internally constructed
	std::shared_ptr<NSPintrct::IntrctMol> imol_ { nullptr }; //*< the simulated system, externally constructed and supplied
	std::shared_ptr<ShakeBonds> shakebds_ { nullptr }; //*< shake bonds determined based on bonds in the simulated moleculare system
	std::shared_ptr<StochasticDynamics::State> state_ { nullptr }; //*<current SD-related states
	std::shared_ptr<StochasticDynamics::State> buffstate_ { nullptr };
	std::shared_ptr<NSPdstl::RandomEngine<>> rng_ { nullptr }; //*<random engine used by SD
	std::shared_ptr<std::vector<std::vector<int>>> temperaturegroups_ { nullptr }; //*< atoms can be divided into groups, each group has its own bath temperature
	std::vector<double> bathtemperatures_;
	std::vector<double> vmasses_;
	std::vector<bool> isfixed_; //<*whether an atom position should be fixed in the SD simulation
	std::vector<std::shared_ptr<ResTraint> > restraints_;
	bool restrainhelices_ { false };
	int nfixedcrd_ { 0 };
	int nstepsrun_ { 0 };
	int nblsteps_ { 50 };
	int sscodesteps_ { 100 }; //<* every number of steps, secondary structure codes is recalculated using the new coordinates.
	bool shakeon_ { true }; //<*whether to constrain bond lengths
	SDRunPara initpara_;
	void init_temperaturegrps(const SDRunPara &para);
};
std::shared_ptr<MolSDRun> newmolsdrun(const SDRunPara &sdpara,
		const NSPintrct::MolSystmPara &mpara,
		const NSPintrct::IntrctPara &ipara);
MolSDRun makemolsdrun(const SDRunPara &sdpara,
		const NSPintrct::MolSystmPara &mpara,
		const NSPintrct::IntrctPara &ipara);
}
#endif /* MOLSDRUN_H_ */

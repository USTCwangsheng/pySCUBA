/*
 * runsddefaults.h
 *
 *  Created on: 2019年12月18日
 *      Author: hyiu
 */
#include "sampling/runsddefaults.h"
#include "geometry/quatfit.h"
#include <memory>
using namespace NSPsampling;
using namespace NSPintrct;
struct RMSDCalculator {
	//double operator()(const std::shared_ptr<std::vector<NSPgeometry::XYZ>> & crd1,
	//		const std::shared_ptr<std::vector<NSPgeometry::XYZ>> & crd2) const {
	double operator()(const SDCallBack_default::Config &cf1,
			const SDCallBack_default::Config &cf2) const {
		return sqrt(NSPgeometry::QuatFit().setup(*(cf1.crd), *(cf2.crd)));
	}
};
SDCallBack_default::SDCallBack_default(const SDRunPara &sdpara) {
	printsteps = sdpara.printsteps;
	savepdbsteps = sdpara.savepdbsteps;
	doannealing = sdpara.doannealing;
	this->pdbfilename = sdpara.outpdbfile;
	if (doannealing) {
		double th = sdpara.annealingscheme[0];
		double tl = sdpara.annealingscheme[1];
		int agrp=sdpara.annealinggroup;
		int nstep_cycle = (int) sdpara.annealingscheme[2];
		int step_h = (int) sdpara.annealingscheme[3];
		int step_d = (int) sdpara.annealingscheme[4];
		tascheme = TemperatureAnnealing(th, tl, nstep_cycle, step_h, step_d,agrp);
	}
}
bool SDCallBack_default::operator()(MolSDRun &run) const {
	std::array<double, NSPintrct::IntrctBlck::ENESIZE> energies;
	double epot_tot = run.imol()->sum_energies(energies);
	if (run.nstepsrun() == 0) {
		epot_av = epot_tot;
		storedconfig.N = run.initpara().storetopn;
		storedconfig.similaritycut = run.initpara().nr_rmsd;
		storedconfig.stored.clear();
	} else {
		double a = run.initpara().epot_decay;
		epot_av = epot_av * a + epot_tot * (1.0 - a);
	}
	if (run.nstepsrun() == 0 || (run.nstepsrun() + 1) % printsteps == 0) {
		double temp = run.temperature();
		std::cout << "Step " << run.nstepsrun() << " Temperature= " << temp
				<< "  Pot_tot= " << epot_tot << " Restraint_ene: ";
		for (auto &r : run.restraints()) {
			std::cout<<r->restriantname() << ":" << r->ene()<<"  ";
		}
		std::cout << std::endl;
	}
	if ((run.nstepsrun() + 1) % savepdbsteps == 0) {
		std::ofstream ofs(pdbfilename);
		run.imol()->writepdb(ofs);
		ofs.close();
	}
	if (doannealing) {
		if(tascheme.annealinggrp<0)
			run.changetemperature(tascheme.temperature(run.nstepsrun() + 1));
		else
			run.changetemperature(tascheme.temperature(run.nstepsrun() + 1),tascheme.annealinggrp);
	}
	if (epot_av < run.initpara().emaxforstore
			&& storedconfig.trystore(epot_av)) {
		auto pcrd = std::shared_ptr<std::vector<NSPgeometry::XYZ>>(
				new std::vector<NSPgeometry::XYZ>(
						run.imol()->recollectcrds_all()));
		Config cf;
		cf.crd = pcrd;
		cf.energies = energies;
		cf.cfstep.push_back(run.nstepsrun());
		cf.pot_tot.push_back(epot_tot);/////
		std::vector<std::string> restrainterms;
		for (auto &r : run.restraints()) {
			cf.restraint_enes.push_back(r->ene());
			restrainterms.push_back(r->restriantname());
			//std::cout << r->restriantname() << std::endl;
		}
		if (storedconfig.store(cf, epot_av, RMSDCalculator())
				&& run.nstepsrun() - lastwritetopstep > 500) {
			lastwritetopstep = run.nstepsrun();
			std::cout << " Top energies: ";
			for (int i = 0; i < storedconfig.nstored(); ++i) {
				std::cout << "  " << storedconfig.getstored(i).score;
				run.imol()->changecrds(*(storedconfig.getstored(i).object.crd));
				std::ofstream ofs(run.initpara().topnpdbname + std::to_string(i) + ".pdb");
				ofs << "COMMENT Pot_tot= ";
				for (auto e : storedconfig.getstored(i).object.pot_tot)
					ofs << e << " ";
				for (auto e : storedconfig.getstored(i).object.cfstep)
				ofs <<"Saved config in step :"<< e <<"  Non-redundant RMSD standard: >"<<run.initpara().nr_rmsd*10<<"  Energy requirement: <"<< run.initpara().emaxforstore << std::endl;
				
				ofs << "COMMENT internal energies= ";
				const std::string str[] = { "BOND","ANGLE","IMPDIH","STERIC","PHIPSI","LOCALSTR","LOCALHB","SITEPAIR","ROTAMER","SCPACKING" };
				int nen = 0;
				for (auto e : storedconfig.getstored(i).object.energies)
					ofs <<str[nen++]<<":"<< e << " ";
				ofs << std::endl;
				int resname = 0;
				ofs << "COMMENT restraint energies= ";				
				for (auto e : storedconfig.getstored(i).object.restraint_enes)
					ofs << restrainterms[resname++] << ":" << e << " ";
				ofs << std::endl;
				run.imol()->writepdb(ofs);
				ofs.close();
			}
			std::cout << std::endl;
			run.imol()->changecrds(*pcrd);
		}
	}
	return false;
}
void NSPsampling::readparameters(const std::string &parfile,
		MolSystmPara &mpara, IntrctPara &ipara, SDRunPara &spara,
		std::string &jobname) {
	NSPdataio::ControlFile cf;
	cf.readfile(parfile);
	mpara = makemolsystmparam(std::string(), cf);
	if (jobname == "auto") {
		jobname = mpara.jobname;
	} else {
		mpara.jobname = jobname;
	}
	ipara = makeintrctparam(jobname, cf);
	spara = makesdrunparam(jobname, cf);
	return;
}


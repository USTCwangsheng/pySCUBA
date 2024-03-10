/*
 * runsddefaults.h
 *
 *  Created on: 2019年12月18日
 *      Author: hyiu
 */

#ifndef RUNSDDEFAULTS_H_
#define RUNSDDEFAULTS_H_
#include "sampling/molsdrun.h"
#include "iblock/intrctmol.h"
#include "iblock/molsystmpara.h"
#include "iblock/intrctparam.h"
#include "sampling/sdrunpara.h"
#include "sampling/temperatureannealing.h"
#include "dstl/nrtopn.h"
namespace NSPsampling{
class SDCallBack_default{
public:
	int printsteps{-1};
	int savepdbsteps{-1};
	bool doannealing{false};
	TemperatureAnnealing tascheme; //*< temperature annealing scheme
	std::string pdbfilename{"sdout.pdb"};
	 struct Config{
			std::shared_ptr<std::vector<NSPgeometry::XYZ>> crd;
			std::array<double,NSPintrct::IntrctBlck::ENESIZE> energies;
			std::vector<double> restraint_enes;
			std::vector<double> pot_tot;///////////添加
			std::vector<int> cfstep;
	};
//	mutable NSPdstl::NRTopN<std::shared_ptr<std::vector<NSPgeometry::XYZ>>> storedconfig;
	mutable NSPdstl::NRTopN<Config> storedconfig;
	mutable double epot_av{0.0};
	mutable int lastwritetopstep{0};
	SDCallBack_default(const SDRunPara &sdpara);
	bool operator() (MolSDRun &run) const;
};


void readparameters(const std::string &parfile,
		NSPintrct::MolSystmPara &mpara,
		NSPintrct::IntrctPara &ipara, SDRunPara &spara,
		std::string &jobname);
}




#endif /* RUNSDDEFAULTS_H_ */

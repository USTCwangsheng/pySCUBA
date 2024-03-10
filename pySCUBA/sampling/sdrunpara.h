/*
 * sdrunpara.h
 *
 *  Created on: 2019年12月16日
 *      Author: hyiu
 */

#ifndef SDRUNPARA_H_
#define SDRUNPARA_H_
#include "iblock/blkselector.h"
#include "dataio/controlfile.h"
#include <vector>
namespace NSPsampling{
struct SDRunPara{
	std::string jobname{""};
	std::string outpdbfile{"sdout.pdb"};
	std::string topnpdbname{ "enetop_" };
	int storetopn{50}; //store topn non-redundant, lowest potential energy configurations
	double emaxforstore{0.0};
	double nr_rmsd{0.2}; //RMSD cutoff for non-redundant configurations
	double epot_decay{0.5};
//	std::string rmsdrestraintfile{""};
	std::string restraintsfile{""};
	bool doshake{true};
	bool shakeinitcrd{true};
//	bool chngintrctbyss{false}; //manipulate the MC interactions according to the intended secd. str. sequence
//	bool avoidlongstrand{false};
//	std::vector<std::string> intendedss;  //intended secondary structure sequence for each chain
	int newnbliststeps{50};
	int newsscodesteps{100};
	int printsteps{100};
	int savepdbsteps{200};
//	int maxstrandlength{100};
//	int ncheckstrandsteps{200};
	double gamma{0.1};
	double atommass{14.};
	double timestep{0.002};
	unsigned randomseed{5};
/**
 * control radius of gyration restraint
 */
//	double krg_rest{0.0}; //<* force constant for radius of gyration restraint
//	double rg0{1.0}; //<* turning point for switching on the radius of gyration restraint
//	std::map<int,std::set<int>> rg_residues; //<* residues of each each chain involved in radius of gyration restraint
/*
 * temperature control
 */
	bool doannealing{false};
	mutable std::vector<std::vector<int>> atomgroups;
	std::vector<std::string> temperaturegrps;
	std::vector<double> temperatures;
	std::vector<double> annealingscheme; //<* parameters specifying temperature annealing schemes
	int annealinggroup{-1}; //<* the temperature group to which annealing will be applied
	SDRunPara(){;}
	SDRunPara(const std::vector<std::string> &controllines);
};
inline SDRunPara makesdrunparam(const std::string &jobname,const NSPdataio::ControlFile &cf,
		const std::string & controlname="SDRunPar"){
		std::vector<std::string> lines=cf.getcontrolines(controlname);
		lines.push_back("JobName = "+jobname);
		return SDRunPara(lines);
}
}



#endif /* SDRUNPARA_H_ */

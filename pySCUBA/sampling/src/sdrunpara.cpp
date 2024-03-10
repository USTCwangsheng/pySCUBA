/*
 * sdrunpara.cpp
 *
 *  Created on: 2019年12月16日
 *      Author: hyiu
 */
#include "sampling/sdrunpara.h"
#include "dataio/parameters.h"
#include "iblock/intrctparam.h"
using namespace NSPsampling;
SDRunPara::SDRunPara(const std::vector<std::string> &controllines){
	   std::map<std::string, double> doublepars {{"GAMMA",0.1},
//	   {"TimeStep",0.002},{"K_RgRestraint",0.0},{"Rg0_Restraint",1.0}};
	   {"TimeStep",0.002}};
		std::map<std::string, std::vector<std::string>> stringvecpars {
//			{"TemperatureGroups",{"ALL"}},{"IntendedSS",{}}
			{"TemperatureGroups",{"ALL"}}
		};
		std::map<std::string, std::vector<double>> doublevecpars {
			{"GroupTemperatures",{1.0}},
			{"AnnealingScheme",{1.0,1.0,1000000,1000000,0}},
			{"StoreTopConfig",{50,2.0,0.5,-1000.0}}
		};
		std::map<std::string, std::string> stringpars {{"JobName",""},
//			{"OutPDBFile","sdout.pdb"},{"Rg_Residues",""},{"RMSDRestraintFile",""}};
		{"OutPDBFile","sdout.pdb"},{"RestraintsFile",""},{"TopnPdbName","enetop_"} };
//		std::map<std::string, std::vector<int>> intvecpars {{"AvoidLongStrand",{0,100,200}}};
		std::map<std::string, std::vector<int>> intvecpars {};
		std::map<std::string, int> intpars{{"PrintParameters",0},
			{"RecalcNeighborListSteps",50},{"RecalcSSSteps",500},
			{"PrintSteps", 100},{"SavePDBSteps",200},
//			{"RandomSeed",5},{"DOShake",1},{"DOAnnealing",0},{"ChngIntrctBySS",0},{"AnnealingGroup",-1}};
			{"RandomSeed",5},{"DOShake",1},{"DOAnnealing",0},{"AnnealingGroup",-1}};
	NSPdataio::ParameterSet pset;
	pset.initdefaultkeys(doublepars,stringpars,intpars,doublevecpars,stringvecpars,intvecpars);
	pset.adjustvalues(controllines);
	pset.getval("GAMMA",&(this->gamma));
	if (gamma < 0.01 || gamma>20) {
		std::cout << "gamma="<< gamma <<".too large or too small" << std::endl;
		exit(0);
	}
	pset.getval("TimeStep",&(this->timestep));
	if (timestep < 0.0005 || timestep>0.02) {
		std::cout << "timestep=" << timestep << ".too large or too small" << std::endl;
		exit(0);
	}
	pset.getval("OutPDBFile",&(this->outpdbfile));
	pset.getval("TopnPdbName", &(this->topnpdbname));
	pset.getval("PrintSteps",&(this->printsteps));
	pset.getval("SavePDBSteps",&(this->savepdbsteps));
	pset.getval("RecalcNeighborListSteps",&(this->newnbliststeps));
	if (newnbliststeps < 5 || newnbliststeps>500) {
		std::cout << "newnbliststeps=" << newnbliststeps << ".too large or too small" << std::endl;
		exit(0);
	}
	pset.getval("RecalcSSSteps",&(this->newsscodesteps));
	if (newsscodesteps < 5 || newsscodesteps>5000) {
		std::cout << "newsscodesteps=" << newsscodesteps << ".too large or too small" << std::endl;
		exit(0);
	}
	std::vector<double> storescheme;
    pset.getval("StoreTopConfig",&storescheme);
    this->storetopn=storescheme[0];
    this->nr_rmsd=storescheme[1]*A2NM;
	if (storescheme[1] < 0.1 || storescheme[1]>10) {
		std::cout << "StoreTopConfig:too large or too small RMSD cutoff for non-redundant configurations" << std::endl;
		exit(0);
	}
    this->epot_decay=storescheme[2];
    this->emaxforstore=storescheme[3];
	/*pset.getval("K_RgRestraint",&(this->krg_rest));
	if(krg_rest>0){
		pset.getval("Rg0_Restraint", &(this->rg0));
		std::string rgresidues;
		pset.getval("Rg_Residues", &rgresidues);
		if(!rgresidues.empty()){
			this->rg_residues=NSPintrct::BlkSelector(rgresidues).selectedblks;
		}
	}*/
//	std::vector<std::string> tmpgrps;
	pset.getval("TemperatureGroups",&(this->temperaturegrps));
	pset.getval("GroupTemperatures",&(this->temperatures));
//	assert(tmpgrps.size()==temperatures.size());
//	if(tmpgrps[0]=="ALL") this->atomgroups.clear();
	for(auto &t:temperatures) t*=KBT;
	pset.getval("JobName",&(this->jobname));
	int seed;
	pset.getval("RandomSeed",&seed);
	this->randomseed=seed;
	int doshake;
	pset.getval("DOShake",&doshake);
	if (doshake != 1 && doshake != 0) { std::cout << "doshake=" << doshake << ".Invalid values" << std::endl;
	exit(0);
	}
	this->doshake=(doshake!=0);
	if(!this->doshake) this->shakeinitcrd=false;
	int annealing;
	pset.getval("DOAnnealing", &annealing);
	if (annealing != 1 && annealing != 0 && annealing!=-1) {
		std::cout << "annealing=" << annealing << ".Invalid values" << std::endl;
		exit(0);
	}
	this->doannealing=(annealing !=0);
	if(doannealing){
		pset.getval("AnnealingScheme", &(this->annealingscheme));
		pset.getval("AnnealingGroup", &(this->annealinggroup));
		if ( this->annealinggroup ==-1 && this->temperaturegrps[0] !="ALL"){
			std::cout <<"Must specify the temperature group to which annealing will be applied"<<std::endl;
			exit(1);
		}
	}

//	std::vector<std::string> iss;
//	pset.getval("IntendedSS", &iss);
//	if(!iss.empty()) this->intendedss=iss;
//	int chngintrct;
//	pset.getval("ChngIntrctBySS",&chngintrct);
//	pset.getval("RMSDRestraintFile", &rmsdrestraintfile);
	pset.getval("RestraintsFile", &restraintsfile);
//	std::cout <<"intendedss =" << iss[0] <<" chngintrtct =" <<chngintrct<<std::endl;
//	if(chngintrct ==1) this->chngintrctbyss=true;
/*	std::vector<int> avoidls;
	pset.getval("AvoidLongStrand",&avoidls);
	if(avoidls[0] !=0) {
		this->avoidlongstrand=true;
		this->maxstrandlength=avoidls[1];
		this->ncheckstrandsteps=avoidls[2];
	}*/
	std::string fixedresidues;
/*	pset.getval("FixedResidues",&fixedresidues);
	NSPintrct::BlkSelector sel(fixedresidues);
	this->allfixedblks=sel.selectedblks;
	std::string mcfixed;
	pset.getval("MainChainFixedResidues",&mcfixed);
	NSPintrct::BlkSelector selmc(mcfixed);
	this->mainchainfixedblks=selmc.selectedblks;*/
	int printpara;
	pset.getval("PrintParameters",&printpara);
	if(printpara !=0){
		auto & oss=NSPintrct::Results::ofstreams(this->jobname,"_par",".dat");
		oss<<"START SDRunPar"<<std::endl;
		pset.readline("PrintParameters = 0");
		pset.printparameters(oss);
		oss<<"END SDRunPar"<<std::endl;
	}
}

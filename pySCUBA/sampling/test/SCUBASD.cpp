/*
 * testsdrun.cpp
 *
 *  Created on: 2019年12月16日
 *      Author: hyiu
 */
#include "sampling/runsddefaults.h"
using namespace NSPintrct;
using namespace NSPsampling;
/**
 * This program illustrates how to use the MolSDRun class to carry out SD simulations
 */
int main(int argc, char **argv){
	std::cout <<"usage: program parameterfile nsteps" <<std::endl;
	std::cout <<"Will run nsteps SD simulations for a system defined by the parameter file"<<std::endl;
	std::cout <<"The simulated structure will be saved regularly to the outpdbfile during the simulation (overwritten mode)"<<std::endl;
	std::cout <<"Temperature, total potential energy at regular intervals: "<<std::endl;
	assert(argc>=3);
	std::string parfile(argv[1]);
	int nsteps=std::stoi(std::string(argv[2]));
	std::string jobname("auto");
	MolSystmPara mpara;
	IntrctPara ipara;
	SDRunPara spara;
	readparameters(parfile,mpara,ipara,spara,jobname);
	SDCallBack_default callback(spara);
	MolSDRun sdrun=makemolsdrun(spara,mpara,ipara);
	sdrun.runsteps(nsteps,callback);
}





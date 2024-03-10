/*
 * mkparfiletemplates.cpp
 *
 *  Created on: 2019年12月12日
 *      Author: hyiu
 */
#include "iblock/intrctparam.h"
#include "sampling/sdrunpara.h"
#include "sampling/loopsampler.h"
#include "iblock/molsystmpara.h"
using namespace NSPintrct;
/**
 * This program generates input file template which can be edited
 * to give proper input files for respective SCUBA programs.
 */
int main(){
	std::vector<std::string> control;
	control.push_back("JobName = tmplt");
	control.push_back("PrintParameters = 1");
    MolSystmPara molpara(control);
	NSPsampling::SDRunPara parasdrun(control);
	IntrctPara param(control);
	NSPsampling::LpSamplerPara lpspara(control);
}





/*
 * printpartemp.cpp
 *
 *  Created on: 2020年8月13日
 *      Author: hyiu
 */
#include "iblock/molsystmpara.h"
#include "iblock/intrctparam.h"
#include "sampling/sdrunpara.h"

using namespace NSPintrct;
using namespace NSPsampling;
/**
 * This program will printout defined parameter keys and values, which may be used as a template  for
 * generate the true parameters
 */
int main(int argc, char **argv){
	std::vector<std::string> lines;
	lines.push_back("PrintParameters = 1");
	MolSystmPara mpara(lines);
	IntrctPara ipara(lines);
	SDRunPara spara(lines);
}


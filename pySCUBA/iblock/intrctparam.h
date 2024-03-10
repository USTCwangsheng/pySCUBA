/*
 * intrctparam.h
 *
 *  Created on: 2019年12月4日
 *      Author: hyiu
 */

#ifndef IBLOCK_INTRCTPARAM_H_
#define IBLOCK_INTRCTPARAM_H_
#include "iblock/blcktopo.h"
#include "dataio/controlfile.h"
#include "iblock/vsctypes.h"
#include <iostream>
#include <fstream>
namespace NSPintrct{

/**
 * An utility to help specifying and open files to store various results
 */
struct Results{
	 static std::ofstream &ofstreams(
			 const std::string &jobname,const std::string &outtype,
			 const std::string &ext=".out");
};

/**
 * Parameters controlling the calculation of forces and energies.
 *
 * An object is usually constructed by reading data from a control file,
 * so that the user can control the calculations by specifying paramters in the control file.
 */
struct IntrctPara{
	IntrctPara(const std::vector<std::string> & intrctcontrolines);
	IntrctPara(){;}
	bool enedetails{true};
	bool phipsicodes_updateonly{false};
	bool sscodes_updateonly{false};
	bool calc_sscodes{true};
	bool calc_nblists{true};
	std::string jobname{"j0"};
	double mm_neighbor_cutoff2 {0.25}; //in nm2;
	double ms_neighbor_cutoff2 {0.64};
	double ss_neighbor_cutoff2 {0.64};
	double pl_neighbor_cutoff2 {0.64};
	double ll_neighbor_cutoff2 {0.64};
	double max_cutoff2{0.64};
	double weight_coval{1.0};
	double weight_steric{1.0};
	double weight_phipsi{1.0};
	double weight_localstr{1.0};
	double weight_sitepair{1.0};
	double weight_rotamer{1.0};
	double weight_scpacking{1.0};
	double weight_localhb{1.0};
};

/**
 * read parameters from a control file
 */
inline IntrctPara makeintrctparam(const std::string &jobname,const NSPdataio::ControlFile &cf,
		const std::string & controlname="InteractionPar"){
		std::vector<std::string> lines=cf.getcontrolines(controlname);
		lines.push_back("JobName = "+jobname);
		return IntrctPara(lines);
}

/**
 * Some energy weights (read from the control file) are rescaled (to change unit from KBT
 * for actual energy unit) before use.
 */
inline IntrctPara mkrescaledparam(const IntrctPara &inipar,double s=KBT){
	IntrctPara res=inipar;
	res.weight_phipsi *=s;
	res.weight_rotamer *=s;
	res.weight_localstr *=s;
	res.weight_sitepair *=s;
	res.weight_localhb *=s;
	return res;
}
}



#endif /* IBLOCK_INTRCTPARAM_H_ */

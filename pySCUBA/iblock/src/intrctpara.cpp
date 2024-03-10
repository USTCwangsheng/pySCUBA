/*
 * intrctpara.cpp
 *
 *  Created on: 2019年12月4日
 *      Author: hyiu
 */
#include "iblock/intrctparam.h"
#include "dataio/parameters.h"
#include <map>
using namespace NSPintrct;
std::ofstream & Results::ofstreams(const std::string &jobname,
		const std::string &outtype, const std::string &ext){
	static std::map<std::string,std::ofstream> streams;
	std::string filename=jobname+outtype+ext;
	if(streams.find(filename)==streams.end()){
		streams[filename]=std::ofstream(filename);
	}
	return streams[filename];
}
IntrctPara::IntrctPara(const std::vector<std::string> & intrctcontrollines){
	   std::map<std::string, double> doublepars { {"CovalentWeight",1.0},
		      { "MCStericWeight", 1.0 }, {
				"PhiPsiWeight", 1.0 }, { "RotamerWeight", 1.0 },
				{ "LocalStrWeight", 1.0 }, { "SitePairWeight", 1.0 }, {"SCPackingWeight",1.0},
				{ "LocalHBWeight", 1.0 }};
		std::map<std::string, std::vector<std::string>> stringvecpars {};
		std::map<std::string, std::vector<double>> doublevecpars {};
		std::map<std::string, std::string> stringpars {{"JobName",""}};
		std::map<std::string, std::vector<int>> intvecpars {};
		std::map<std::string, int> intpars{{"PrintParameters",0},{"WriteEneDetails",0}};
	NSPdataio::ParameterSet pset;
	pset.initdefaultkeys(doublepars,stringpars,intpars,doublevecpars,stringvecpars,intvecpars);
	pset.adjustvalues(intrctcontrollines);
	pset.getval("CovalentWeight",&(this->weight_coval));
	pset.getval("MCStericWeight",&(this->weight_steric));
	pset.getval("PhiPsiWeight",&(this->weight_phipsi));
	pset.getval("LocalStrWeight",&(this->weight_localstr));
	pset.getval("RotamerWeight",&(this->weight_rotamer));
	pset.getval("SitePairWeight",&(this->weight_sitepair));
	pset.getval("SCPackingWeight",&(this->weight_scpacking));
	pset.getval("LocalHBWeight",&(this->weight_localhb));
	pset.getval("JobName",&(this->jobname));
	int edetails;
	pset.getval("WriteEneDetails",&edetails);
	this->enedetails=(edetails!=0);
	int printpara;
	pset.getval("PrintParameters",&printpara);
	if(printpara !=0){
		auto & oss=Results::ofstreams(this->jobname,"_par",".dat");
		oss<<"START InteractionPar"<<std::endl;
		pset.readline("PrintParameters = 0");
		pset.printparameters(oss);
		oss<<"END InteractionPar"<<std::endl;
	}
}



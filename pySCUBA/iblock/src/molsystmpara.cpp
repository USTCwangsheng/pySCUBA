/*
 * molsystmpara.cpp
 *
 *  Created on: 2019年12月17日
 *      Author: hyiu
 */
#include "iblock/intrctparam.h"
#include "iblock/molsystmpara.h"
#include "dataio/parameters.h"
#include "iblock/intrctmol.h"
//#include "fullsite/fullsite.h"
using namespace NSPintrct;

MolSystmPara::MolSystmPara(const std::vector<std::string> & controllines){
	   std::map<std::string, double> doublepars {};
		std::map<std::string, std::vector<std::string>> stringvecpars {
			{"Sequence1L",std::vector<std::string>()}};
		std::map<std::string, std::vector<double>> doublevecpars {};
		std::map<std::string, std::string> stringpars {{"JobName",""},{"PDBStart",""},
			{"MainChainFixedResidues",""},{"FixedResidues",""},
			{"ActiveResidues",""},{"SideChainActiveResidues",""},
			{"SoftSideChainResidues",""}};
		std::map<std::string, std::vector<int>> intvecpars {};
		std::map<std::string, int> intpars{{"PrintParameters",0},{"ChangePDBSequence",0}};
	NSPdataio::ParameterSet pset;
	pset.initdefaultkeys(doublepars,stringpars,intpars,doublevecpars,stringvecpars,intvecpars);
	pset.adjustvalues(controllines);
	pset.getval("JobName",&(this->jobname));
	pset.getval("Sequence1L",&(this->sequence1l));
	pset.getval("PDBStart",&(this->pdbstart));
	std::string mcf;
	pset.getval("MainChainFixedResidues",&mcf);
	if(!mcf.empty()) this->fixedmainchainresidues=BlkSelector(mcf).selectedblks;
	std::string fixed;
	pset.getval("FixedResidues",&fixed);
	if(!fixed.empty()) this->fixedresidues=BlkSelector(fixed).selectedblks;
	std::string soft;
	pset.getval("SoftSideChainResidues",&soft);
	if(!soft.empty()) this->softsidechainresidues=BlkSelector(soft).selectedblks;
	std::string allactive;
	pset.getval("ActiveResidues",&allactive);
	if(!allactive.empty()) this->activeresidues=BlkSelector(allactive).selectedblks;
	std::string scactive;
	pset.getval("SideChainActiveResidues",&scactive);
	if(!scactive.empty()) this->sidechainactiveresidues=BlkSelector(scactive).selectedblks;
	if(!fixedresidues.empty()||!fixedmainchainresidues.empty()){
		assert(activeresidues.empty() && sidechainactiveresidues.empty());
	}
	int change;
	pset.getval("ChangePDBSequence",&change);
	if(change!=0) this->changepdbsequence=true;
	if(changepdbsequence){
		assert(!sequence1l.empty() && !pdbstart.empty());
	}
	if(jobname.empty())jobname=pdbstart;
	int printpara;
	pset.getval("PrintParameters",&printpara);
	if(printpara !=0){
		auto & oss=Results::ofstreams(this->jobname,"_par",".dat");
		oss<<"START MolSystmPar"<<std::endl;
		pset.readline("PrintParameters = 0");
		pset.printparameters(oss);
		oss<<"END MolSystmPar"<<std::endl;
	}
}
using namespace NSPproteinrep;
std::shared_ptr<IntrctMol> NSPintrct::make_molsystm(const MolSystmPara &molpara,
		const IntrctPara &ipara){
	std::string pdbfilename=molpara.pdbstart;
/*	std::vector<std::vector<FullSite>> chains=
			readfullsitesfrompdb(pdbfilename,true);*/
	NSPproteinrep::AAConformersInModel model;
	model.readpdbfile(pdbfilename);
	auto & chains=model.conformers;
	std::shared_ptr<IntrctMol> imol=std::shared_ptr<IntrctMol>(new IntrctMol(model));
	for(int c=0;c<chains.size();++c){
		for(int r=0;r<chains[c].size()-1;++r){
			double omiga=chains[c][r].omega_c(chains[c][r+1]);
			if(omiga>-1.57 && omiga <1.57) {
				std::cout <<"Found cis peptide bond after residue " <<c<<"-"<<r<<std::endl;
				(*imol)(NSPdstl::Idx2D{c,r}).precis=true;
			}
		}
	}
	for(auto &csel:molpara.softsidechainresidues){
			for(auto &r:csel.second){
				(*imol)(NSPdstl::Idx2D{csel.first,r}).softsc=true;
			}
	}
	bool fixspecified=!molpara.fixedmainchainresidues.empty() || ! molpara.fixedresidues.empty();
	bool activespecified=!molpara.activeresidues.empty() || !molpara.sidechainactiveresidues.empty();
	if(fixspecified &&activespecified){
		std::cout <<"Simultaneous specification of fixed and active residue selections in MolSystmPara not supported"<<std::endl;
		exit(1);
	}
	if (fixspecified){
		imol->specifyfixedblks(molpara.fixedresidues,molpara.fixedmainchainresidues);
	}
	if(activespecified ){
		imol->specifyactiveblks(molpara.activeresidues,molpara.sidechainactiveresidues);
	}
	imol->my_intrctpara=std::shared_ptr<IntrctPara>(new IntrctPara(mkrescaledparam(ipara)));
	return imol;
}

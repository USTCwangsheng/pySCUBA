/*
 * ligpepsampler.cpp
 *
 *  Created on: 2020年10月23日
 *      Author: hyiu
 */
#include "sampling/ligpepsampler.h"
using namespace NSPsampling;
using namespace NSPintrct;
#include "dataio/parameters.h"
LigPSamplerPara::LigPSamplerPara(const std::vector<std::string> &controllines){
	   std::map<std::string, double> doublepars {{"ENEAverageDecay",0.99},
	   {"ENEVarStop",10.0}};
		std::map<std::string, std::vector<std::string>> stringvecpars {{"InterfaceResidues",{}}};
		std::map<std::string, std::vector<double>> doublevecpars {{"ContactRestraint",{0.0,0.0}}};
		std::map<std::string, std::string> stringpars {{"JobName",""},
		{"PDBStart",""},{"SaveConfigFile",""}};
		std::map<std::string, std::vector<int>> intvecpars {};
		std::map<std::string, int> intpars{{"MaxSDSteps",0},{"RandomSeed",1357}
			,{"PrintParameters",0},{"Verbose",1},{"LigPepLength",0}};
	NSPdataio::ParameterSet pset;
	pset.initdefaultkeys(doublepars,stringpars,intpars,doublevecpars,stringvecpars,intvecpars);
	pset.adjustvalues(controllines);
	pset.getval("MaxSDSteps",&(this->maxsdsteps));
	pset.getval("ENEAverageDecay",&(this->enedecay));
	pset.getval("ENEVarStop",&(this->sdenevarcut));
	std::vector<double> contactres;
	pset.getval("ContactRestraint",&contactres);
	this->ncontacts0=contactres[0];
	this->kcontactres=contactres[1];
	pset.getval("JobName",&(this->jobname));
	pset.getval("PDBStart",&(this->pdbstart));
	pset.getval("SaveConfigFile", &(this->saveconfigfile));
	pset.getval("RandomSeed", &(this->randomseed));
	pset.getval("LigPepLength",&(this->ligpeplength));
	std::vector<std::string> ifselect;
		pset.getval("InterfaceResidues",&ifselect);
	for(auto & s:ifselect){
		this->ifresidues = ifresidues +" " + s;
	}
	int printpara;
	pset.getval("PrintParameters",&printpara);
	if(printpara !=0){
		auto & oss=NSPintrct::Results::ofstreams(this->jobname,"_par",".dat");
		oss<<"START LigPepSamplerPar"<<std::endl;
		pset.readline("PrintParameters = 0");
		pset.printparameters(oss);
		oss<<"END LigPepSamplerPar"<<std::endl;
	}
	int verbosein;
	pset.getval("Verbose",&verbosein);
	verbose=true;
	if(verbosein==0) verbose=false;
}
void NSPsampling::readpara_ligpepsampler(const std::string &parfile,
		LigPSamplerPara &lpspara, NSPintrct::IntrctPara &ipara, SDRunPara &spara,
		std::string &jobname){
	NSPdataio::ControlFile cf;
		cf.readfile(parfile);
		lpspara=makeligpsamplerparam(std::string(),cf);
		if(jobname=="auto"){
			jobname=lpspara.jobname;
		} else {
			lpspara.jobname=jobname;
		}
		ipara=makeintrctparam(jobname,cf);
		spara=makesdrunparam(jobname,cf);
		return;
}
LigPepSampler NSPsampling::mkligpepsampler(const LigPSamplerPara &lpspara,
		const NSPintrct::IntrctPara &ipara,
		const SDRunPara &sdpara){
	NSPproteinrep::AAConformersInModel model;
	model.readpdbfile(lpspara.pdbstart);
	IntrctMol mol0(model);
	LigPepSampler res;
	res.mypara()=lpspara;
	res.mypara().sdprintsteps=sdpara.printsteps;
	res.setupmol(mol0);
	res.imol()->my_intrctpara=std::shared_ptr<IntrctPara>(new IntrctPara(mkrescaledparam(ipara)));
	res.setupsdrun(sdpara);
	return res;
}
void LigPepSampler::setupmol(const IntrctMol &receptor){
	imol_=std::shared_ptr<IntrctMol>(new IntrctMol(receptor));
	BlkSelector ifresidues(mypara_.ifresidues);
	 std::vector<NSPdstl::Idx2D>  interfaceresidues;
	 for(auto &sc:ifresidues.selectedblks){
		 for(auto &r:sc.second)
			 interfaceresidues.push_back({sc.first,r});
	 }
	 NSPdstl::RandomEngine<>::getinstance().reseed(mypara_.randomseed);
	 this->ligandpep_=std::shared_ptr<LigandPep> (new LigandPep(*imol_,interfaceresidues));
	this->ligandpep_->addligandpep(imol_.get(),mypara_.ligpeplength,false);
}
void  LigPepSampler::setupsdrun(const SDRunPara &sdpara){
	   BlkSelector ligblks;
	   double ctgdmin{ 0.5 };
	   double ctgdsmall{ 0.8 };
	   double ctgdoff{ 2.0 };
	   int ligcid=imol_->nchains()-1; //chain id of ligand peptide
	   ligblks.selectedblks.insert(std::make_pair(ligcid,std::set<int>()));
	   for(int i=0;i<imol_->nresidues(ligcid);++i){
		   ligblks.selectedblks[ligcid].insert(i);
	   }
		imol_->specifyactiveblks(ligblks.selectedblks,typename NSPintrct::BlkSelector::SelectedBlks());
		optimizer_.init(sdpara,imol_,false);
		BlkSelector ifresidues(mypara_.ifresidues);
		std::shared_ptr<ResTraint> cr=std::shared_ptr<ResTraint>(new
				ContactRestraint(*imol_,ifresidues,ligblks,mypara_.ncontacts0,mypara_.kcontactres,ctgdmin, ctgdsmall, ctgdoff));
		optimizer_.addrestraint(cr);
	}
void LigPepSampler::buildandoptimize(){
	   this->rebuildligand();
		const std::vector<double> &initcrd=imol_->getcrds_all_d();
		optimizer_.reinitstate(initcrd,0);
		SDCallBackOptSD callbck(mypara_);
		optimizer_.runsteps(mypara_.maxsdsteps,callbck);
}
void LigPepSampler::saveconfig(std::ofstream &ofs) const {
	std::array<double,IntrctBlck::ENESIZE>  energies;
	double etot_=imol_->sum_energies(energies);
	ofs<<etot_<<std::endl;
	for(auto e:energies) ofs <<" "<<e;
	ofs<<std::endl;
	for(auto &r:optimizer_.restraints()) ofs<<r->ene();
	ofs<<std::endl;
	const std::vector<NSPgeometry::XYZ> & xyz= imol_->recollectcrds_all();
	int ligcid=imol_->nchains()-1;
	int begin=imol_->getblck({ligcid,0}).aoffset;
	int nligatoms=xyz.size()-begin;
	assert(nligatoms=mypara_.ligpeplength*4);
	ofs<<nligatoms<<std::endl;
	for(int i=begin;i<xyz.size();++i){
		ofs<<xyz.at(i).toString()<<std::endl;
	}
}

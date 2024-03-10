/*
 * molsdrun.cpp
 *
 *  Created on: 2019年12月16日
 *      Author: hyiu
 */
#include "sampling/molsdrun.h"
#include "iblock/intrctmol.h"
#include "sampling/rmsdrestraintreader.h"
#include <algorithm>
#include <cfloat>
#include <memory>
using namespace NSPsampling;
void MolSDRun::init_temperaturegrps(const SDRunPara &para){
	this->temperaturegroups_=std::shared_ptr<std::vector<std::vector<int>>>(new std::vector<std::vector<int>>());
	if(para.temperaturegrps[0]=="ALL") return;
	if(para.temperaturegrps[0]=="MC+SC"){
		std::vector<int> mcatoms=NSPintrct::BlkSelector().selectmcatoms(*(this->imol_));
		std::vector<int> scatoms=NSPintrct::BlkSelector().selectscatoms(*(this->imol_));
		this->temperaturegroups_->push_back(mcatoms);
		this->temperaturegroups_->push_back(scatoms);
	}
	if(para.temperaturegrps[0]=="SELECTIONS"){
		//int n = (para.temperaturegrps.size() + 1) / 2;
		std::string tpgline{""};
		for (int c = 1;c < para.temperaturegrps.size();++c) {
			//std::cout << para.temperaturegrps[c] << std::endl;
			tpgline= tpgline + para.temperaturegrps[c] +" ";
		}

		std::vector<std::string> temp_word;
		NSPutils::split(tpgline,temp_word,";");

		int n = temp_word.size();
		for(int i=1;i<= n;++i){
			std::string temperaturegrpn = temp_word[i-1];
			//std::cout << i << temperaturegrpn << std::endl;
			NSPintrct::BlkSelector selector(temperaturegrpn);
			std::cout << "temperaturegrps:  "<<temperaturegrpn <<" in temperature :"<<(para.temperatures[i-1])/2.4942 <<"*KBT"<< std::endl;
			std::vector<int> atoms=selector.selectatoms(*(this->imol_));
			this->temperaturegroups_->push_back(atoms);
		}
	}
	//check if every atom has a temperature
	std::vector<int> atomlabels(this->imol_->natoms(),0);
	for(auto &g:*(this->temperaturegroups_)){
		for(auto i:g) atomlabels[i]=1;
	}
	int nsum=0;
	for(int i:atomlabels) nsum+=i;
	if(nsum != this->imol_->natoms()){
		std::cout<<"Not all atoms have been assigned into a temperature group. Check TemperatureGroups parameter"<<std::endl;
		exit(1);
	}
	return;
	std::cout<<"Unknown specification of TemperatureGroups "<<para.temperaturegrps[0] <<std::endl;
	exit(1);
}
bool MolSDRun::init(const SDRunPara &para,
		std::shared_ptr<NSPintrct::IntrctMol> imol,bool shakeinit) {
	imol_ = imol;
	 initpara_=para;
	std::vector<NSPgeometry::XYZ> initcrds = imol->recollectcrds_all();
	// calculate potential energies at starting point
/*	std::vector<NSPgeometry::XYZ> dedx;
	imol_->forces_all(initcrds, dedx);
	imol_->my_intrctpara.phipsicodes_updateonly=true;
	imol_->my_intrctpara.sscodes_updateonly=true;
	imol_->my_intrctpara.calc_nblists=false;
	imol_->my_intrctpara.calc_sscodes=false;*/
	isfixed_ = imol_->isfixedvec();
	initrandomengine<unsigned int>(para.randomseed);
	sd_ = std::shared_ptr<StochasticDynamics>(new StochasticDynamics());
	this->init_temperaturegrps(initpara_);
	para.atomgroups=*(this->temperaturegroups_);
	sd_->init(para, *imol);
/*	if(para.krg_rest>0){
              std::vector<double> wrg;
              std::vector<int> rgatoms=imol_->atomindices(para.rg_residues);
              if(!rgatoms.empty()){
            	  wrg.resize(imol_->natoms(),0.0);
            	  for(auto a:rgatoms) wrg[a]=1.0;
              }
              std::shared_ptr<ResTraint> rgrest=std::shared_ptr<ResTraint>(
            		  new RgRestraint(para.rg0*A2NM,para.krg_rest*KBT,wrg));
              this->addrestraint(rgrest);
	}*/
//	RMSDRestraintReader().addstructrestraints(*this);
	if (!para.restraintsfile.empty()) {
		restraints_ = readrestraints(*this->imol_, para.restraintsfile);
		std::cout << "num of restraints:" << restraints_.size() << std::endl;
	}
//	chngintrctbyss();
	//print restraint information
	if(!restraints_.empty()){
		std::cout <<"The following restraints will be applied: "<<std::endl;
		for(auto &r:restraints_){
               r->printinfo(std::cout);
		}
	}
	this->bathtemperatures_=para.temperatures;
//	this->temperaturegroups_=std::shared_ptr<std::vector<std::vector<int>>>(new std::vector<std::vector<int>>());
//	for(auto &l:para.atomgroups) temperaturegroups_->push_back(l);
	shakebds_ = make_shakebonds(*imol);
	shakeon_ = para.doshake;
//	restrainhelices_=para.restrainhelices;
	shakebds_->seton(shakeon_);
	nblsteps_ = para.newnbliststeps;
	sscodesteps_ = para.newsscodesteps;
	vmasses_ = sd_->masses();
	int natoms = imol_->natoms();
	nfixedcrd_ = 0;
	for (int i = 0; i < natoms; ++i) {
		if (isfixed_[i]) {
			for (int idx = 3 * i; idx < 3 * i + 3; ++idx) {
				vmasses_[idx] = MASS_MAX;
				++nfixedcrd_;
			}
		}
	}

	state_ = sd_->make_initstate(NSPgeometry::XYZvtodoublev(initcrds), *rng_,
			&vmasses_);
	buffstate_ = std::shared_ptr<StochasticDynamics::State>(
			new StochasticDynamics::State);
	*buffstate_ = *state_;
	nstepsrun_ = 0;
	if(!shakeinit) return true;
	return sd_->shakestate(*state_, shakebds_->shakefunction(),
			para.shakeinitcrd, &vmasses_);

}
/*void MolSDRun::chngintrctbyss(){
    if (!this->initpara_.chngintrctbyss) return;
	auto & ssseqs=this->initpara_.intendedss;
//	std::cout <<"intendedss =" << initpara_.intendedss[0] <<" chngintrtct =" <<initpara_.chngintrctbyss<<std::endl;
	if(ssseqs.empty()) return;
	std::shared_ptr<HelixRestraint>  hrest=std::shared_ptr<HelixRestraint>(new HelixRestraint(this->imol_.get()));
	int cid=0;
	for(auto &ssvec:ssseqs){
		std::cout <<"chain "<<cid <<" IntendedSS:" <<ssvec<<std::endl;
		for(int i=0;i<ssvec.size();++i){
			auto &blk1= this->imol_->getblck(NSPdstl::Idx2D(cid,i));
			 if(ssvec[i]=='C'){
				blk1.weight_mc=0.1;
			 } else if (ssvec[i]=='H'){
				 if(i+4<ssvec.size() ){
					 int a1=blk1.aoffset+blk1.natoms()-1;
					 int a2=this->imol_->getblck(NSPdstl::Idx2D(cid,i+4)).aoffset;
					 if(ssvec[i+4]=='H') hrest->hbrestraints.push_back(DisRestraint(a1,a2,0.36,0.4,5000));
				 }
			 }
		}
		++cid;
	}
	this->restraints_.push_back(hrest);
}*/
std::shared_ptr<MolSDRun> NSPsampling::newmolsdrun(const SDRunPara &sdpara,
		const NSPintrct::MolSystmPara &mpara, const NSPintrct::IntrctPara &ipara){
	std::shared_ptr<MolSDRun> sdrun=std::shared_ptr<MolSDRun>(new MolSDRun());
	std::shared_ptr<NSPintrct::IntrctMol> imol=NSPintrct::make_molsystm(mpara,ipara);
	sdrun->init(sdpara,imol);
	return sdrun;
}
MolSDRun NSPsampling::makemolsdrun(const SDRunPara &sdpara,
		const NSPintrct::MolSystmPara &mpara, const NSPintrct::IntrctPara &ipara){
	MolSDRun sdrun;
	std::shared_ptr<NSPintrct::IntrctMol> imol=NSPintrct::make_molsystm(mpara,ipara);
	sdrun.init(sdpara,imol);
	return sdrun;
}
/*
void MolSDRun::chngstrandflankingweights(){
	static const double flankingw=0.1;
	int nchains= imol_->nchains();
	for(int cid=0;cid<nchains; ++cid){
		for (int i=0;i<imol_->sscodes[cid].size();++i){
			if(!this->initpara_.intendedss.empty()){
				if(this->initpara_.intendedss[cid][i]=='C') continue;
			}
			 imol_->getblck(NSPdstl::Idx2D(cid,i)).weight_mc=1.0;
		}
		std::vector<std::pair<int,int>> helices,strands;
		imol_->ssregions(cid,helices,strands);
		for(auto &s:strands){
	//		std::cout <<"strand length " <<s.second<<std::endl;
			if(s.second>=this->initpara_.maxstrandlength){
				 std::cout <<"Long strand in chain " <<cid <<" start at "<< s.first <<" length " <<s.second<<std::endl;
				 for(int i=1;i<3;++i){
					  std::array<int,2> idxs;
					  idxs[0]=s.first-i;
					  idxs[1]=s.first+s.second+i-1;
					  for (auto idx:idxs ){
						       if(idx<0 || idx >=imol_->sscodes.at(cid).size()) continue;
						       if(!this->initpara_.intendedss.empty()){
						    	   if(this->initpara_.intendedss[cid][idx]=='C') continue;
						       }
						     imol_->getblck(NSPdstl::Idx2D(cid,idx)).weight_mc=flankingw;
					  }
				 }
			}
		}
	}
}
*/

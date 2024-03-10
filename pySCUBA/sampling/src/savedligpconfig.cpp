/*
 * savedligpconfig.cpp
 *
 *  Created on: 2020年10月25日
 *      Author: hyiu
 */
#include "sampling/savedligpconfig.h"

using namespace NSPsampling;
using namespace NSPintrct;
using namespace NSPproteinrep;
void SavedLigPConfig::init(LigPSamplerPara &lpars) {
	std::cout <<"Starting setup "<<std::endl;
	AAConformersInModel model;
	model.readpdbfile(lpars.pdbstart);
	IntrctMol mol0(model);
	LigPepSampler res;
	res.mypara() = lpars;
	res.setupmol(mol0);
	this->imol = res.imol();
	this->peplength = lpars.ligpeplength;
	savedconfigfile = lpars.saveconfigfile;
	std::cout <<"Finished setup "<<std::endl;
	this->readallconfigs();
}

SavedLigPConfig::Config SavedLigPConfig::readnextconfig(std::istream &is) {
	int configid;
	char buffer[120];
//	std::cout <<"Start reading Config " <<std::endl;
	try {
		is.getline(buffer, 120);
//		std::cout <<"firstline "<< buffer<< std::endl;
		if (!is.good())
			return SavedLigPConfig::Config();
		std::vector<std::string> words = NSPdataio::parseline(
				std::string(buffer), std::vector<int>());
		SavedLigPConfig::Config res;
		res.id = std::stoi(words[0]);
		res.epot_tot = std::stod(words[1]);
//		std::cout <<res.id<<" "<<res.epot_tot<<std::endl;
		is.getline(buffer, 120);
		words = NSPdataio::parseline(std::string(buffer), std::vector<int>());
		for (int i = 0; i < NSPintrct::IntrctBlck::ENESIZE; i++) {
			try {
				res.energies[i] = std::stod(words[i]);
			} catch (...) {
				res.energies[i] = 0.0;
			}
		}
		is.getline(buffer, 120);
		words = NSPdataio::parseline(std::string(buffer), std::vector<int>());
		res.e_contactrest = std::stod(words[0]);
		is.getline(buffer, 120);
		words = NSPdataio::parseline(std::string(buffer), std::vector<int>());
		int natoms = std::stoi(words[0]);
		assert(natoms = 4 * (this->peplength));
//		std::cout <<res.e_contactrest<<" "<<natoms<<std::endl;
		res.crds.resize(natoms);
		for (int i = 0; i < natoms; ++i) {
			is.getline(buffer, 120);
			words = NSPdataio::parseline(std::string(buffer),
					std::vector<int>());
			res.crds[i].x_ = std::stod(words[0]);
			res.crds[i].y_ = std::stod(words[1]);
			res.crds[i].z_ = std::stod(words[2]);
		}
//		std::cout <<"Read Config " <<res.id<<std::endl;
		return res;
	} catch (...) {
		return Config();
	}
}
void SavedLigPConfig::readallconfigs() {
	std::ifstream ifs(this->savedconfigfile);
	configs.clear();
	Config next = readnextconfig(ifs);
	while (!next.crds.empty()) {
		configs.push_back(next);
		next = readnextconfig(ifs);
	}
}
void SavedLigPConfig::determinetopn(int n, double rmsdcut){
		CompareConfig cc;
		cc.slp=this;
		topnconfigs=NSPdstl::NRTopN<int>(n,rmsdcut);
		for(int i=0;i<configs.size();++i){
		    double et=configs[i].epot_tot+configs[i].e_contactrest;
			if(topnconfigs.trystore(et)){
				 topnconfigs.store(configs[i].id,et,cc);
			}
		}
}
void SavedLigPConfig::writepdb_topn(int n,double rmsdcut){
	   if(topnconfigs.nstored()<n) determinetopn(n,rmsdcut);
	   int idoffset=configs[0].id;
	   for(int i=0;i<n;++i){
		   if(i>=topnconfigs.nstored()) break;
		   auto x=topnconfigs.getstored(i);
		   int id=x.object;
		   Config & cf=configs[id-idoffset];
		   std::string pdbfile="top_"+std::to_string(i)+".pdb";
		   std::ofstream ofs(pdbfile);
		   writepdb(ofs,cf);
		   ofs.close();
	   }
}
void  SavedLigPConfig::writepdb(std::ofstream &ofs, const Config &cf,bool writeene){
	std::vector<NSPgeometry::XYZ> xyz=imol->recollectcrds_all();
	int aoffset=imol->getblck({imol->nchains()-1,0}).aoffset;
	std::copy(cf.crds.begin(),cf.crds.end(),xyz.begin()+aoffset);
	imol->changecrds(xyz);
	if(writeene){
		ofs <<"E_tot = " << cf.epot_tot <<", E_contactrestraint=  "<<cf.e_contactrest<<std::endl;
	}
	imol->writepdb(ofs);
}

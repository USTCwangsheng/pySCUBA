/*
 * savedloops.cpp
 *
 *  Created on: 2019年12月25日
 *      Author: hyiu
 */

#include "sampling/savedloops.h"
#include "dstl/sortindex.h"
#include "dataio/inputlines.h"
#include "geometry/calculators.h"
using namespace NSPsampling;
using namespace NSPintrct;
SavedLoops::Config SavedLoops::readnextconfig(std::istream &is) {
	int configid;
	char buffer[120];
	try {
		is.getline(buffer,120);
//		std::cout <<"firstline "<< buffer<< std::endl;
		if(!is.good()) return Config();
		std::vector<std::string> words=NSPdataio::parseline(std::string(buffer),std::vector<int>());
		Config res;
		res.id=std::stoi(words[0]);
		is.getline(buffer,120);
		words=NSPdataio::parseline(std::string(buffer),std::vector<int>());
		for (int i = 0; i < NSPintrct::IntrctBlck::ENESIZE; i++) {
			   	try{
			   			res.energies[i]=std::stod(words[i]);
			   	} catch(...)
			   	{
			   		res.energies[i]=0.0;
			   	}
//			   		std::cout <<words[i]<<std::endl;
		}
//		std::cout <<"num loop atoms: " <<nloopatoms<<std::endl;
		res.crds.resize(nloopatoms);
		for (int i = 0; i < nloopatoms; ++i) {
			is.getline(buffer,120);
			words=NSPdataio::parseline(std::string(buffer),std::vector<int>());
			res.crds[i].x_ =std::stod(words[0]);
			res.crds[i].y_ =std::stod(words[1]);
			res.crds[i].z_=std::stod(words[2]);
//			std::cout <<i<< " "<<res.crds[i].toString()<<std::endl;
		}
		return res;
	} catch (...) {
		return Config();
	}
}
void SavedLoops::readallconfigs() {
	std::ifstream ifs(this->loopfile);
	configs.clear();
	Config next = readnextconfig(ifs);
	while (!next.crds.empty()) {
		configs.push_back(next);
		next = readnextconfig(ifs);
	}
}
void SavedLoops::init(LpSamplerPara &para) {
	NSPproteinrep::AAConformersInModel model;
	model.readpdbfile(para.pdbstart);
//	std::vector<std::vector<NSPproteinrep::FullSite>> pdbchains =
//			NSPproteinrep::readfullsitesfrompdb(para.pdbstart);
	IntrctMol mol0(model);
	LoopSampler res;
	res.mypara()=para;
	for(auto & fp:para.flankingpositions){
		int cid=mol0.mappdbkeyint()->chainNumber(fp.first.first);
		int left=mol0.mappdbkeyint()->
				posiNumber(NSPproteinrep::PdbReader::reskey(fp.first.second),fp.first.first);
		int right=mol0.mappdbkeyint()->
				posiNumber(NSPproteinrep::PdbReader::reskey(fp.second.second),fp.first.first);
		res.mypara().flankingsites.push_back({{cid,left},{cid,right}});
	}
	res.setupmol(mol0, res.mypara().flankingsites,res.mypara().looplengths);
	imol=res.imol();
	nloopatoms = 0;
	nloops=para.nloops;
	looplengths=para.looplengths;
    reindexedflankingsites=res.reindexedflankingsites();
    loopfile=para.outloopfile;
	for (int l=0;l<nloops;l++){
		for(int i=0;i<looplengths[l];i++){
				int cid=reindexedflankingsites[l].first.first;
				int rid=reindexedflankingsites[l].first.second+1+i;
				nloopatoms+=(*imol)({cid,rid}).natoms();
		}
	}
	copypdbconfig();
	std::cout <<"Copied loop config from pdb read in" <<std::endl;
	this->readstartconfig();
	std::cout<<"RMSD between startconfig and pdbconfig: "<<std::endl;
	std::cout <<rmsd(pdbconfig,startconfig)/A2NM <<std::endl;
	this->readallconfigs();
}
void SavedLoops::writermsdene(const std::string &filename){
	std::vector<std::pair<double,double>> rmsdene;
	std::vector<int> ids;
	for(auto & c:configs){
		double et=0.0;
		for(auto &e:c.energies) et+=e;
		if(et>100.0) continue;
	    double d=rmsd(c,pdbconfig)/A2NM;
	     rmsdene.push_back({d,et});
	     ids.push_back(c.id);
	}
	NSPdstl::sort12(rmsdene,ids);
	std::ofstream ofs(filename);
	for(int i=0;i<rmsdene.size();++i){
		ofs <<ids[i] <<"\t" << rmsdene[i].first <<"\t"<<rmsdene[i].second<<std::endl;
	}
	ofs.close();
	ofs.open("lowestrmsd.pdb");
	writepdb(ofs,configs[ids[0]-configs[0].id],true);
}
void SavedLoops::writepdb(std::ofstream & ofs,const Config & cf,bool writeene){
	int cidx=0;
	auto &crds=cf.crds;
    if(writeene){
    	double etot=0.0;
    	for(auto &e:cf.energies) etot+=e;
    	ofs<<"REMARK Etot =" <<etot<<std::endl<<"REMARK  Ecomponents:";
//    	for(auto &e:cf.energies) ofs <<" " <<e;
    	std::vector<std::string> EneName {"BOND", "ANGLE", "IMPDIH", "STERIC", "PHIPSI",
    	"LOCALSTR", "LOCALHB", "SITEPAIR", "ROTAMER", "SCPACKING"};
    	int k = 0;
    	for(auto &e:cf.energies)
    	{
    		ofs << " " << EneName[k] << ":" << e;
    		k++;
    	}
    	ofs<<std::endl;
    	double d=rmsd(cf,pdbconfig)/A2NM;
    	ofs<<"REMARK RMSD from native: " <<d <<std::endl;
    }
	for(int l=0;l<looplengths.size();l++){
		for(int i=0;i<looplengths[l];i++){
			int cid=reindexedflankingsites[l].first.first;
			int rid=reindexedflankingsites[l].first.second+1+i;
			int nat=(*imol)({cid,rid}).natoms();
			(*imol)({cid,rid}).crds.resize(nat);
			std::copy(crds.begin()+cidx,crds.begin()+cidx+nat,(*imol)({cid,rid}).crds.begin());
			cidx+=nat;
		}
	}
	(*imol).writepdb(ofs);
}
SavedLoops::Config SavedLoops:: copypdbconfig(){
	pdbconfig.crds.clear();
	for(int l=0;l<looplengths.size();l++){
			for(int i=0;i<looplengths[l];i++){
				int cid=reindexedflankingsites[l].first.first;
				int rid=reindexedflankingsites[l].first.second+1+i;
				int nat=(*imol)({cid,rid}).natoms();
				for(int a=0;a<nat;++a){
					pdbconfig.crds.push_back((*imol)({cid,rid}).crds.at(a));
				}
			}
	}
	return pdbconfig;
}
double SavedLoops::rmsd(const Config & config1, const Config & config2) {
	double d2=0.0;
	int natoms=config1.crds.size();
	for(int a=0;a<natoms;++a){
		d2 += (config1.crds[a]-config2.crds[a]).squarednorm();
	}
	return sqrt(d2/(double) natoms);
}
double SavedLoops::pprmsd(const Config & config1, const Config & config2, std::vector<int> seq) {
	double d2=0.0;
	double d1 = 0.0; // sum of pp;
	int nres=config1.crds.size()/4;

	for(int a=0;a<nres;++a){
		if (a > 0)
		{
			double phi1 = torsion(config1.crds[(a-1)*4+seq[2]],config1.crds[a*4+seq[0]],config1.crds[a*4+seq[1]],config1.crds[a*4+seq[2]]);
			double phi2 = torsion(config2.crds[(a-1)*4+seq[2]],config2.crds[a*4+seq[0]],config2.crds[a*4+seq[1]],config2.crds[a*4+seq[2]]);
			double dphi = 360*(phi1-phi2)/3.14159265;
			if (dphi > 360) dphi -= 360;
			else if (dphi < -360) dphi += 360;
			d2 += dphi*dphi;
			d1 += abs(dphi);
		}
		if (a < nres-1)
		{
			double psi1 = torsion(config1.crds[a*4+seq[0]],config1.crds[a*4+seq[1]],config1.crds[a*4+seq[2]],config1.crds[(a+1)*4+seq[0]]);
			double psi2 = torsion(config2.crds[a*4+seq[0]],config2.crds[a*4+seq[1]],config2.crds[a*4+seq[2]],config2.crds[(a+1)*4+seq[0]]);
			double dpsi = 360*(psi1-psi2)/3.14159265;
			if (dpsi > 360) dpsi -= 360;
			else if (dpsi < -360) dpsi += 360;
			d2 += dpsi*dpsi;
			d1 += abs(dpsi);
		}
	}

	return sqrt(d2/(double) (2*nres-2)); // each res have 1 phi & 1 psi.
}

void SavedLoops::determinetopn(int n, double rmsdcut){
		CompareConfig cc;
		cc.slp=this;
		topnconfigs=NSPdstl::NRTopN<int>(n,rmsdcut);
		for(int i=0;i<configs.size();++i){
			double et=0.0;
			for(auto &e:configs[i].energies) et +=e;
			if(topnconfigs.trystore(et)){
				 topnconfigs.store(configs[i].id,et,cc);
			}
		}
}
void SavedLoops::determinetopnPP(int n, double rmsdcut){
		CompareConfigPP cc;
		cc.slp=this;
		topnconfigs=NSPdstl::NRTopN<int>(n,rmsdcut);
		for(int i=0;i<configs.size();++i){
			double et=0.0;
			for(auto &e:configs[i].energies) et +=e;
			if(topnconfigs.trystore(et)){
				 topnconfigs.store(configs[i].id,et,cc);
			}
		}
}
void SavedLoops::writepdb_topn(int n,double rmsdcut){
	   if(topnconfigs.nstored()<n) determinetopn(n,rmsdcut);
	   int idoffset=configs[0].id;
	   for(int i=0;i<n;++i){
		   if(i>=topnconfigs.nstored()) break;
		   auto x=topnconfigs.getstored(i);
		   int id=x.object;
		   double et=x.score;
		   std::string pdbfile="top_"+std::to_string(i)+".pdb";
		   std::ofstream ofs(pdbfile);
		   writepdb(ofs,configs[id-idoffset]);
		   ofs.close();
	   }
}
void SavedLoops::savenonredundant(const std::string &file,
		const std::string & rmsdenefile,double rmsdcut){
	topnconfigs.stored.clear();
	this->determinetopn(configs.size(),rmsdcut);
//	this->determinetopnPP(configs.size(),rmsdcut);
	std::cout<<"Total number of non redundant configs: " <<topnconfigs.nstored()<<std::endl;
	std::ofstream ofs(file);
	std::ofstream ofs2(rmsdenefile);
	for(int i=0;i<topnconfigs.nstored();++i){
		auto x=topnconfigs.getstored(i);
		Config &c=configs[x.object-configs[0].id];
		double d=rmsd(c,pdbconfig);
		ofs<<i  <<" " <<c.id<<" " <<x.score<<" " <<d <<std::endl;
		for (auto e:c.energies)
			ofs <<" " <<e;
		ofs <<std::endl;
		for (auto &x:c.crds) ofs <<x.toString() <<std::endl;
		ofs2<<d/A2NM<<"\t" <<x.score<<std::endl;
	}
}

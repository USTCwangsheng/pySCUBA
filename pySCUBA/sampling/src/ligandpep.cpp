/*
 * ligandpep.cpp
 *
 *  Created on: 2020年10月17日
 *      Author: hyiu
 */
#include "sampling/ligandpep.h"
#include "dstl/sortindex.h"
#include "dstl/randomengine.h"
#include "backbone/backbonebuilder.h"
#include "iblock/molmodeler.h"
using namespace NSPgeometry;
using namespace NSPsampling;
using namespace NSPintrct;
using namespace NSPproteinrep;
#define INITDISTANCE 1.0  //distance in nanometer
void LigandPep::setup(const IntrctMol & receptor,
	            const std::vector<NSPdstl::Idx2D> & interfaceresidues){
	      const std::vector<XYZ> & rcrds= receptor.recollectcrds_all();
	      this->receptorcenter_=XYZ(0,0,0);
	      for(auto &a:rcrds){
	    	  this->receptorcenter_ =this->receptorcenter_+a;
	      }
	      this->receptorcenter_=(1.0/(double) rcrds.size())*this->receptorcenter_;
	      for(auto & idx:interfaceresidues){
	    	  interfacecacrds_.push_back(receptor.getblck(idx).getcrd("CA"));
	      }
	      ifrespairs_.clear();
	      if(interfacecacrds_.size()<2) return;
	      std::vector<double> rij_ca;
	      std::vector<std::pair<int,int>> ij_idx;
	      for(int i=0;i<interfacecacrds_.size()-1;++i){
	    	  for(int j=i+1; j<interfacecacrds_.size();++j){
	    		  rij_ca.push_back(-(interfacecacrds_[i]-interfacecacrds_[j]).squarednorm());
	    		  ij_idx.push_back({i,j});
	    	  }
	      }
	      NSPdstl::sort12(rij_ca,ij_idx);
	      int idxm=ij_idx.size();
	      if(idxm>3) idxm=3+(idxm-3)/3;
	      ifrespairs_.resize(idxm);
	      std::copy(ij_idx.begin(),ij_idx.begin()+idxm,ifrespairs_.begin());
}
XYZ LigandPep::randomligandcenter(
			XYZ *direction) const{
	XYZ result;
	auto &rng=NSPdstl::RandomEngine<>::getinstance().realrng(0.0,1.0);
	auto &ing=NSPdstl::RandomEngine<>::getinstance().intrng(0,ifrespairs_.size());
	if(ifrespairs_.empty()){
		XYZ rd=(interfacecacrds_[0]-this->receptorcenter_);
		rd=(1.0/sqrt(rd.squarednorm()))*rd;
		result=interfacecacrds_[0]+ INITDISTANCE*rd + XYZ(rng,0.1);
		XYZ rr(rng,1.0);
		*direction=cross(rd,rr);
		return result;
	} else {
		int idx=ing();
		XYZ rm=0.5*(interfacecacrds_[ifrespairs_[idx].first]+interfacecacrds_[ifrespairs_[idx].second]);
		XYZ rd=rm-this->receptorcenter_;
		rd =(1.0/sqrt(rd.squarednorm()))*rd;
		result=rm+INITDISTANCE*rd+XYZ(rng,0.1);
		if(rng()<0.5) *direction=interfacecacrds_[ifrespairs_[idx].second]-interfacecacrds_[ifrespairs_[idx].first];
		else  *direction=interfacecacrds_[ifrespairs_[idx].first]-interfacecacrds_[ifrespairs_[idx].second];
		*direction=*direction+XYZ(rng,0.1*(sqrt(direction->squarednorm())));
		return result;
	}
}
std::vector<NSPgeometry::XYZ> LigandPep::genligandcrds(int ligandlength) const{
	XYZ dir;
	auto center=(1.0/A2NM)*this->randomligandcenter(&dir);
	BackBoneSite bs0;
	genbackbonesite(nullptr,false,-120.0,120.0,&bs0);
	std::vector<NSPproteinrep::BackBoneSite> bss=
				BackBoneBuilder::buildforwardbackbone(ligandlength,
				bs0,std::vector<std::pair<int,int>>(),
				std::vector<std::pair<int,int>> (), std::set<int>() );
		BackBoneBuilder::movechainto(XYZ(0.0,0.0,0.0),dir,true,bss);
		std::vector<NSPgeometry::XYZ> xyz;
		for(auto & bs:bss){
			for(int i=0;i<4;++i) xyz.push_back(bs.getcrd(BackBoneSite::NCRD+3*i));
		}
		XYZ oldc(0.0,0.0,0.0);
		for(auto & r:xyz) oldc=oldc+r;
		oldc=center-(1/(4.0*ligandlength))*oldc;
		for(auto &r:xyz) r=A2NM*(r+oldc);
		return xyz;
}
void LigandPep::addligandpep(NSPintrct::IntrctMol *imol,int ligandlength,bool gencrd)  const {
	MolModeler modeler;
	modeler.settarget(*imol);
	std::vector<IntrctBlck> blcks;
	for(int i=0;i<ligandlength;++i){
		blcks.push_back(modeler.mkfloatingblck(
				std::make_pair(-10,-10),
				"GLY",BackBoneSite()));
	}
	imol->addchain(blcks);  //the blcks do not have correct coordindates
	if(gencrd) {
		int cid =imol->nchains()-1;
		this->genligandpepcrd(imol,cid);
	}
}
void LigandPep::genligandpepcrd(IntrctMol *imol, int ligandchainid) const {
	int ligandlength=imol->nresidues(ligandchainid);
	std::vector<XYZ> xyz=this->genligandcrds(ligandlength);
	std::vector<XYZ> crd=imol->recollectcrds_all();
	int offset=imol->getblck({ligandchainid,0}).aoffset;
	std::copy(xyz.begin(),xyz.end(),crd.begin()+offset);
	imol->changecrds(crd);
}


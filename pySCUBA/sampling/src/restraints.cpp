/*
 * restraints.cpp
 *
 *  Created on: 2020年8月13日
 *      Author: hyiu
 */
#include "sampling/restraints.h"
#include "geometry/quatfit.h"
using namespace NSPsampling;
using namespace NSPgeometry;

StructRestraint::StructRestraint(const std::vector<XYZ> &crd,
		const std::vector<double> &weights, double kres,double rmsdref,double rmsdswitch,int mode) :
				crdref_(crd),
		weights_(weights), kres_(kres),mode_(mode),rmsdref_(rmsdref) , rmsdswitch_(rmsdswitch){
	assert(crd.size()  == weights_.size());
	wtot_=0.0;
	for(auto w:weights_) wtot_+=w;
}

StructRestraint::StructRestraint(const std::vector<XYZ> &crd,
		const std::vector<std::pair<int,int>> &poslengths, double kres,int mode,double rmsdref) :
		crdref_(crd),	 kres_(kres),mode_(mode),rmsdref_(rmsdref) {
	weights_.assign(crd.size(),0.0);
	wtot_=0.0;
	for(auto iter=poslengths.begin(); iter != poslengths.end(); ++iter){
		for(int i=iter->first; i<iter->first+iter->second;++i){
			weights_[i]=1.0;
			wtot_+=1.0;
		}
	}
}
double StructRestraint::posi_energy(const std::vector<NSPgeometry::XYZ> &crd,std::vector<NSPgeometry::XYZ> *forces) const{
	NSPgeometry::QuatFit qfit;
	std::vector<NSPgeometry::XYZ> crdtmp=crdref_;
	double rmsd2=qfit.fitting(crd,crdtmp,weights_);//changing the refcrd!!!
	const double kresh=0.5*kres_;
	double etot=0.0;
//	double wtot=0.0;
//	double d2=0.0;
	for(int i=0;i<crd.size();++i) {
		if(weights_[i]<1.e-10) continue;
		std::vector<NSPgeometry::XYZ> deriv;
		double b = NSPgeometry::distance(crd[i], crdtmp[i], &deriv);
		if(b<1.e-10) continue;
//		d2+=b*b;
//		wtot+=weights_[i];
		etot += weights_[i]*kresh * b*b;
		double kbdb = -weights_[i]*kres_ * b;
//		if(kbdb >1.e-3 || kbdb<-1.e-3){
//			std::cout <<kbdb<<" "<<kres_<<" "<<b<<" "<<weights_[i]<<std::endl;
//			std::cout<<i <<" " <<crdref_.size()<<std::endl;
//			std::cout<<crd.at(i).toString() <<crdtmp.at(i).toString()<<std::endl;
		forces->at(i) = forces->at(i) + kbdb * deriv[0];
	}
//	double  rmsd2_recalc=d2/wtot;
//	std::cout <<"rmsd2: "<<rmsd2_<< " "<<rmsd2_recalc<<std::endl;
	rmsd2_=rmsd2;
	return etot;
}
double StructRestraint::totalrmsd_energy(const std::vector<NSPgeometry::XYZ> &crd,std::vector<NSPgeometry::XYZ> *forces) const{
	NSPgeometry::QuatFit qfit;
	std::vector<NSPgeometry::XYZ> crdtmp=crdref_;
	double rmsd2=qfit.fitting(crd,crdtmp,weights_);//changing the refcrd!!!
	if(rmsd2<=rmsdref_*rmsdref_) {
		rmsd2_=rmsd2;
		return 0.0;
	}
	double rmsd=sqrt(rmsd2);
	double drmsd=rmsd-rmsdref_;
	double dedrmsd=kres_*drmsd;
	double ene=0.5*dedrmsd*drmsd;
	if(rmsd >rmsdswitch_){
		double d=(rmsdswitch_-rmsdref_);
		dedrmsd=kres_*d;
		ene=0.5*dedrmsd*d+ dedrmsd*(rmsd-rmsdswitch_);
	}
	dedrmsd /=rmsd;
	for(int i=0;i<crd.size();++i) {
		if(weights_[i]<1.e-10) continue;
		forces->at(i) = forces->at(i) - (dedrmsd *weights_[i]/wtot_)* (crd[i]-crdtmp[i]);
	}
//	double  rmsd2_recalc=d2/wtot;
//	std::cout <<"rmsd2: "<<rmsd2_<< " "<<rmsd2_recalc<<std::endl;
	rmsd2_=rmsd2;
	return ene;
}
/*
double RgRestraint::energy(const std::vector<NSPgeometry::XYZ> &crd,
		std::vector<NSPgeometry::XYZ> *forces,const std::vector<double> &w) const {
	if(kres_>0) return energy1(crd,forces,w);
	else return energy2(crd,forces,w);
}
double RgRestraint::energy1(const std::vector<NSPgeometry::XYZ> &crd,
		std::vector<NSPgeometry::XYZ> *forces,const std::vector<double> &w) const{
		if(kres_==0.0) return 0.0;
		NSPgeometry::XYZ center=NSPgeometry::center(crd,w);
		double rtot2=0.0;
		std::vector<NSPgeometry::XYZ> drt2dx;
		std::vector<NSPgeometry::XYZ> deriv;
		int idx=0;
		double wtot=0.0;
		for(auto &c:crd){
			NSPgeometry::XYZ dc=c-center;
			double r2=dc.squarednorm();
			double wgt=1.0;
			if(!w.empty()){
				wgt=w[idx];
			}
			rtot2+=wgt*r2;
			drt2dx.push_back(2.0*wgt*dc);  //drtot2/dx
			wtot +=wgt;
			++idx;
		}
		double rg=sqrt(rtot2/wtot);
		if(rg<=rgbound_) return 0.0;
		double drg=rg-rgbound_;
		double ene=0.5*kres_*drg*drg;
		double dedrtot2=-kres_*drg/(2.0*rg*wtot); //kres*drg * (1/(2*rg)*1/crd.size()
		for(int i=0;i<crd.size();++i){
			(*forces)[i] = (*forces)[i] +dedrtot2*drt2dx[i];
		}
		return ene;
}*/

/*
 * ene=kres_*log(rg^2)
 * kres_: should be proportional to temperature and the number of degrees of freedom
 */
double RgRestraint::logrestr_energy(const std::vector<NSPgeometry::XYZ> &crd,
		std::vector<NSPgeometry::XYZ> *forces) const{
		if(kres_==0.0) return 0.0;
		NSPgeometry::XYZ center=NSPgeometry::center(crd,w_);
		double rtot2=0.0;
		std::vector<NSPgeometry::XYZ> drt2dx;
		std::vector<NSPgeometry::XYZ> deriv;
		int idx=0;
		double wtot=0.0;
		for(auto &c:crd){
			double wgt=1.0;
			if(!w_.empty()){
				wgt=w_[idx];
			}
			NSPgeometry::XYZ dc=c-center;
			double r2=dc.squarednorm();
			rtot2+=wgt*r2;
			wtot +=wgt;
			drt2dx.push_back(2.0*wgt*dc);  //drtot2/dx
			++idx;
		}
		double rg2=rtot2/wtot;
//		std::cout <<"Rg= " << sqrt(rg2) <<std::endl;
		if(rg2<=rgbound_*rgbound_) return 0.0;
		double ene=kres_*log(rg2/(rgbound_*rgbound_));
		double dedrtot2=-kres_/(rg2*wtot); //kres*drg * (1/(2*rg)*1/crd.size()
		for(int i=0;i<crd.size();++i){
			(*forces)[i] = (*forces)[i] +dedrtot2*drt2dx[i];
		}
		return ene;
}
double DisRestraint::energy_attract(const std::vector<NSPgeometry::XYZ> & crd, std::vector<NSPgeometry::XYZ> *forces) const{
	std::vector<NSPgeometry::XYZ> drdx;
	double r=NSPgeometry::distance(crd[a1],crd[a2],&drdx);
//	std::cout <<"r, r0 "<< r <<" "<<r0<<std::endl;
	if (r<=r0) return 0.0;
	double dr=r-r0;
	double ene,dedr;
	if(r<r1) {
		dedr=kres*dr;
		ene=0.5*dedr*dr;
	} else{
		dedr=kres*(r1-r0);
		ene=0.5*dedr*(r1-r0)+dedr*(r-r1);
	}
	(*forces)[a1]=(*forces)[a1] - dedr*(drdx[0]);
	(*forces)[a2]=(*forces)[a2] - dedr*(drdx[1]);
	return ene;
}
double DisRestraint::energy_repulsion(const std::vector<NSPgeometry::XYZ> & crd, std::vector<NSPgeometry::XYZ> *forces) const{
	std::vector<NSPgeometry::XYZ> drdx;
	double r=NSPgeometry::distance(crd[a1],crd[a2],&drdx);
	if (r>=-r0) return 0.0;
	double dr=-r0-r;
	double ene,dedr;
	if( dr>r1){
		dedr=kres*dr;
		ene=0.5*dedr*dr;
	}
	else {
		dedr=kres*(-r0-r1);
		ene=0.5*dedr*(-r0-r1)+dedr*(r1-r);
	}
	(*forces)[a1]=(*forces)[a1] + dedr*(drdx[0]);
	(*forces)[a2]=(*forces)[a2] + dedr*(drdx[1]);
	return ene;
}
double GrpDisRestraint::energy(const std::vector<NSPgeometry::XYZ> &crd, std::vector<NSPgeometry::XYZ> *forces) const{
		ene_ = 0;
		std::vector<NSPgeometry::XYZ> grpcrd(2,{0,0,0});
		std::vector<NSPgeometry::XYZ> grpforces(2,{0,0,0});
		double f1=1.0/(double) grp1.size();
		double f2=1.0/(double) grp2.size();
		for(auto i:grp1) grpcrd[0] = grpcrd[0]+crd.at(i);
		grpcrd[0]=f1*grpcrd[0];
		for(auto i:grp2) grpcrd[1] = grpcrd[1]+crd.at(i);
		grpcrd[1]=f2*grpcrd[1];
		ene_=disres.energy(grpcrd,&grpforces);
		for(auto i:grp1){
			(*forces)[i] = (*forces)[i]+ f1*grpforces[0];
		}
		for(auto i:grp2){
			(*forces)[i]=(*forces)[i] +f2*grpforces[1];
		}
		return ene_;
	}
ContactRestraint::ContactRestraint(const NSPintrct::IntrctMol & imol,
		const NSPintrct::BlkSelector & iresidues_receptor,  const NSPintrct::BlkSelector & iresidues_ligand,
		double nc0_i,double kres_i,double resgdmin_i, double resgdsmall_i, double resgdoff_i): nc0(nc0_i),kres(kres_i), resgdmin(resgdmin_i), resgdsmall(resgdsmall_i), resgdoff(resgdoff_i){
	std::vector<int> ica_receptor;
	std::vector<int> ica_ligand;
	NSPintrct::AtomIdx3D aid{0,0,1};
	for(auto & cb:iresidues_receptor.selectedblks){
		for(auto &r:cb.second){
			ica_receptor.push_back(imol.getblck({cb.first,r}).getatomindex_abs(aid));
		}
	}
	for(auto & cb:iresidues_ligand.selectedblks){
				for(auto &r:cb.second){
					ica_ligand.push_back(imol.getblck({cb.first,r}).getatomindex_abs(aid));
			}
	}
	//add groupcontact between each ligand interfaceresidue CA and all receptor interface residue CAs
	for(auto ica: ica_ligand){
		this->gcs.push_back(GroupContact());
		gcs.back().gdmin = resgdmin_i;
		gcs.back().gdsmall=resgdsmall_i;
		gcs.back().gdoff = resgdoff_i;
		gcs.back().grpd.grps[0]=ica_receptor;
		gcs.back().grpd.grps[1]=std::vector<int>(1,ica);
	}
}
double ContactRestraint::energy(const std::vector<NSPgeometry::XYZ> &crd, std::vector<NSPgeometry::XYZ> *forces) const{
	ene_=0.0;
    double nc=0.0;
    std::vector<NSPintrct::DvDxi> dncdx;
    for(auto &gc:gcs){
    	 nc += gc.theta(crd,&dncdx);
    }
    if(nc < nc0) {
    	ene_=0.5*kres*(nc-nc0)*(nc-nc0);
    	double dednc=kres*(nc-nc0);
    	for(auto &d:dncdx){
    		(*forces)[d.first] =(*forces)[d.first] - d.second*dednc;
    	}
    }
    return ene_;
}
/*HelixRestraint::HelixRestraint(NSPintrct::IntrctMol *im, int flankingmcoff,double kres, double ptarget):
imol(im){
	auto & sscodes=imol->sscodes;
	std::vector<std::string> ssseqs;
	for(auto &cvec:sscodes){
		ssseqs.push_back("");
		for(auto &code:cvec) ssseqs.back().push_back(code.charcode());
	}
	for(auto &ssseq:ssseqs){
		for(int i=1;i<ssseq.size()-1; ++i){
			if(ssseq[i]!='H' && ssseq[i-1]=='H' && ssseq[i+1]=='H') ssseq[i]='H';
		}
	}
	std::vector<std::vector<std::pair<int,int>>> helices;
	for(auto & ssvec:ssseqs){
		int helixstart=-1;
		int helixend =-1;
		int pos=0;
		helices.push_back(std::vector<std::pair<int,int>>());
		for (auto &ss:ssvec){
			if(ss=='H'){
				if(helixstart<0) {
					helixstart=pos;
					helixend=-1;
				}
			} else {
				if(helixstart >=0){
				     helixend=pos-1;
				     helices.back().push_back(std::make_pair(helixstart,helixend));
				     helixstart=-1;
				}
			}
			++pos;
		}
	}
	int cid=0;
	for(auto &hvec:helices){
		for(auto &helix:hvec){
			int len=helix.second-helix.first+1;
			if(len<5){
				for(int i=helix.first; i<=helix.second; ++i) ssseqs[cid][i]='C';
			}
		}
		cid++;
	}
	cid=0;
	for (auto & hvec:helices){
		for(auto &helix:hvec){
			int len=helix.second-helix.first+1;
			if(len<5) continue;
			std::cout <<helix.first<<":"<<helix.second<<std::endl;
//			restraints.push_back(NSPintrct::SSRestraint(cid,helix.first,helix.second,0,ptarget,kres));
			for(int i=helix.first; i<helix.second-2;++i){
				auto & blk1=im->getblck(NSPdstl::Idx2D(cid,i));
				auto &blk2=im->getblck(NSPdstl::Idx2D(cid,i+4));
				int a1=blk1.aoffset+blk1.natoms()-1;
				int a2=blk2.aoffset;
				hbrestraints.push_back(DisRestraint(a1,a2,0.36,0.4,5000));
			}
			if(flankingmcoff>0){
				for(int i=1;i<=flankingmcoff;i++){
					int idx=helix.first-i;
					if(idx>=0) {
						if(ssseqs[cid][idx]!='H') im->getblck(NSPdstl::Idx2D(cid,idx)).weight_mc=0.0;
					}
					idx=helix.second+i;
					if(idx<sscodes[cid].size()){
						if(ssseqs[cid][idx] !='H')
						im->getblck(NSPdstl::Idx2D(cid,idx)).weight_mc=0.0;
					}
				}
			}
		}
		++cid;
		}
//debug
	for(int idx=0;idx<sscodes[0].size();++idx){
	    if(im->getblck(NSPdstl::Idx2D(0,idx)).weight_mc==0.0) std::cout <<"-"<<idx;
	}
	std::cout <<std::endl;
	}
*/

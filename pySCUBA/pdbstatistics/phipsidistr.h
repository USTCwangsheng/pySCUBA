/*
 * phipsidistr.h
 *
 *  Created on: 2016年11月26日
 *      Author: hyliu
 */

#ifndef PDBSTATISTICS_PHIPSIDISTR_H_
#define PDBSTATISTICS_PHIPSIDISTR_H_
#include "dataio/datapaths.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>
namespace NSPpdbstatistics {

class PhiPsiDistr {
public:
	static double coilstep;
	static double glystep;
	static double preprostep;
	static double cisprostep;
	static double transprostep;
	PhiPsiDistr(double step): step_(step){
		phidimension_=(unsigned int) (360.001/step);
		psidimension_=phidimension_;
		distr_.resize(phidimension_*psidimension_,0.005/(double) (phidimension_*psidimension_));
//		distr_.resize(phidimension_*psidimension_,1.e-20);
	}

	static const PhiPsiDistr & coildistr(const std::string & file="coilphipsi.dat") {
		static PhiPsiDistr distr(coilstep);
		if(!distr.distrread_) {
#ifdef _OPENMP
#pragma omp critical(phipsicoil)
			{
				if(!distr.distrread_) distr.readdistr(file);
#pragma omp flush
			}
#else
			distr.readdistr(file);
#endif
		}
		return distr;
	}

	static const PhiPsiDistr & mixcoildistr(double pgly=0.1054,double ppro=0.0623){
		static PhiPsiDistr distr(coilstep);
		if(!distr.distrread_){
#ifdef _OPENMP
#pragma omp critical(phipsimixcoil)
			{
				if(!distr.distrread_) {
#endif
			distr.distrread_=true;
			auto & coil=coildistr();
			auto & gly=glydistr();
			auto & pro=transprodistr();
			auto & prepro=preprodistr();
			distr.distrmax_=-1.0;
			for (int idxphi=0;idxphi<distr.phidimension_;++idxphi){
				for(int idxpsi=0;idxpsi<distr.psidimension_;++idxpsi){
					auto phipsi=distr.getphipsi(idxphi,idxpsi);
					double d=pgly*gly.distr(phipsi.first,phipsi.second) +
							ppro*pro.distr(phipsi.first,phipsi.second) +
							(1-pgly-ppro)*(ppro*prepro.distr(phipsi.first,phipsi.second)+
									(1-ppro)*coil.distr(phipsi.first,phipsi.second));
							d /=(double)(distr.phidimension_*distr.psidimension_);
							if(d>distr.distrmax_) distr.distrmax_=d;
							distr.distr_[idxpsi*distr.phidimension_+idxphi] =d;
				}
			}
#ifdef _OPENMP
#pragma omp flush
			}  //omp distrread
		} //omp critical
#endif
		}
		return distr;
	}

	static const PhiPsiDistr & helixdistr(const std::string & file="helixphipsi.dat") {
		static PhiPsiDistr distr(coilstep);
		if(!distr.distrread_) {
#ifdef _OPENMP
#pragma omp critical(phipsihelix)
			{
				if(!distr.distrread_) distr.readdistr(file);
#pragma omp flush
			}
#else
			distr.readdistr(file);
#endif
		}
		return distr;
	}
	static const PhiPsiDistr & stranddistr(const std::string & file="strandphipsi.dat") {
		static PhiPsiDistr distr(coilstep);
		if(!distr.distrread_) {
#ifdef _OPENMP
#pragma omp critical(phipsistrand)
			{
				if(!distr.distrread_) distr.readdistr(file);
#pragma omp flush
			}
#else
			distr.readdistr(file);
#endif
		}
		return distr;
	}
	static const PhiPsiDistr & glydistr(const std::string &file="glyphipsi.dat") {
		static PhiPsiDistr distr(glystep);
		if(!distr.distrread_) {
#ifdef _OPENMP
#pragma omp critical(phipsigly)
			{
				if(!distr.distrread_) distr.readdistr(file);
#pragma omp flush
			}
#else
			distr.readdistr(file);
#endif
		}
		return distr;
	}
	static const PhiPsiDistr & preprodistr(const std::string & file="preprophipsi.dat") {
		static PhiPsiDistr distr(preprostep);
		if(!distr.distrread_) {
#ifdef _OPENMP
#pragma omp critical(phipsiprepro)
			{
				if(!distr.distrread_) distr.readdistr(file);
#pragma omp flush
			}
#else
			distr.readdistr(file);
#endif
		}
		return distr;
	}
	static const PhiPsiDistr & transprodistr(const std::string & file="transprophipsi.dat") {
		static PhiPsiDistr distr(transprostep);
		if(!distr.distrread_) {
#ifdef _OPENMP
#pragma omp critical(phipsitranspro)
			{
				if(!distr.distrread_) distr.readdistr(file);
#pragma omp flush
			}
#else
			distr.readdistr(file);
#endif
		}
		return distr;
	}
	static const PhiPsiDistr & cisprodistr(const std::string & file="cisprophipsi.dat") {
		static PhiPsiDistr distr(cisprostep);
		if(!distr.distrread_) {
#ifdef _OPENMP
#pragma omp critical(phipsicispro)
			{
				if(!distr.distrread_) distr.readdistr(file);
#pragma omp flush
			}
#else
			distr.readdistr(file);
#endif
		}
		return distr;
	}
	static const PhiPsiDistr &phipsidistr(std::string res, std::string nextres="ALA"){
		if(res=="GLY") return glydistr();
		else if(res=="PRO" || res=="TRANSPRO") return transprodistr();
		else if(res=="CISPRO") return cisprodistr();
		else if(nextres=="PRO") return preprodistr();
		else return coildistr();
	}
	double distr(std::pair<unsigned int,unsigned int> phipsi_indices) const {
		unsigned int idx=phipsi_indices.second*phidimension_+phipsi_indices.first;
		return distr_[idx];
	}
	double distr(double phi,double psi) const {
//		return distr(getindex(phi,psi))*(double)(phidimension_*psidimension_);
		std::pair<unsigned int,unsigned int> bottomleft;
		std::pair<unsigned int,unsigned int> topright;
		std::pair<double,double> xy=getcrdingrid(phi,psi,&bottomleft,&topright);
		double x=xy.first;
		double y=xy.second;
		std::pair<unsigned int,unsigned int> topleft(bottomleft.first,topright.second);
		std::pair<unsigned int,unsigned int> bottomright(topright.first,bottomleft.second);
		double mix=distr(bottomleft)*(1.-x)*(1.-y)+
				   distr(bottomright)*x*(1.-y) +
				   distr(topleft)*(1.-x)*y+
				   distr(topright)*x*y;
		return mix*(double)(phidimension_*psidimension_);
	}
	double mdistr_phi(double phi) const {
//		unsigned int phi_index= getindex(phi,0).first;
		unsigned int left=getindex(phi-0.5*step_,0).first;
		unsigned int right=getindex(phi+0.5*step_,0).first;
		double dphi=phi-getphipsi(left,0).first;
		while(dphi>180.0) dphi-=360.0;
		while(dphi<-180.0) dphi+=360.0;
		double x=dphi/step_;
		double res=0.0;
		for (int i=0;i<psidimension_;++i) {
			res += distr(std::pair<unsigned int,unsigned int>(left,i))*(1.0-x)+
					distr(std::pair<unsigned int,unsigned int>(right,i))*x;
		}
		return res*(double) phidimension_;
	}
	double mdistr_psi(double psi) const {
//		unsigned int psi_index= getindex(0,psi).second;
		unsigned int bottom=getindex(0,psi-0.5*step_).second;
		unsigned int top=getindex(0,psi+0.5*step_).second;
		double dpsi=psi-getphipsi(0,bottom).second;
		while(dpsi>180.0) dpsi-=360.0;
		while(dpsi<-180.0) dpsi+=360.0;
		double y=dpsi/step_;
		double res=0.0;
		for (int i=0;i<phidimension_;++i) {
			res += distr(std::pair<unsigned int,unsigned int>(i,bottom))*(1.-y)+
					distr(std::pair<unsigned int,unsigned int>(i,top))*y;
		}
		return res*(double) psidimension_;
	}
	void estimatedistr(const std::vector<std::pair<double,double>> &phipsi_vals) {
		std::vector<double> counts(phidimension_*psidimension_,0);
		for(auto v:phipsi_vals){
			std::pair<unsigned int,unsigned int> indices=getindex(v.first, v.second);
			counts[indices.second*phidimension_+indices.first] +=1.0;
		}
		distrmax_=-1.0;
		for (unsigned int i=0; i<counts.size();++i) {
			distr_[i]= counts[i]/(double) (phipsi_vals.size());
			if(distr_[i] > distrmax_) distrmax_=distr_[i];
		}
	}
	void writedistr(const std::string &filename,double dcut=1.e-8) const {
		std::ofstream ofs;
		ofs.open(filename.c_str());
		writedistr(ofs,dcut);
		ofs.close();
	}
	void writedistr(std::ostream & os,double dcut=1.e-8) const {
		os <<step_<<std::endl;
		for (unsigned int i=0; i<psidimension_;++i)
			for (unsigned int j=0;j<phidimension_;++j) {
				double d=distr(std::make_pair(j,i));
				if(d < dcut) continue;
				std::pair<double,double> phipsi=getphipsi(j,i);
				os <<  phipsi.first <<"  "<<phipsi.second <<" " <<d <<std::endl;
			}
	}
	void readdistr(const std::string &filename) {
		std::ifstream ifs;
//		std::string path=NSPdataio::datapath()+filename;
		std::string path=NSPdataio::datafilename(filename);
		ifs.open(path.c_str());
		assert(ifs.good());
		readdistr(ifs);
		ifs.close();
	}
	void readdistr(std::istream & is){
		double phi;
		double psi;
		double d;
		distrmax_=-1.0;
		double newstep;
		is >>newstep;
		double mindistr=0.005/(double) (phidimension_*psidimension_);
		if(newstep !=step_) {
			step_=newstep;
			phidimension_=(unsigned int) (360.001/step_);
			psidimension_=phidimension_;
			distr_.resize(phidimension_*psidimension_,mindistr);
		}
//		double scale=1.0/(1+1.e-4);
		while (is.good())  {
			is >> phi >>psi >>d;
			if(is.good()) {
				std::pair<unsigned int, unsigned int> indices=getindex(phi,psi);
				if(d<mindistr) d=mindistr;
				distr_[indices.second*phidimension_+indices.first] = d;
				if(d >distrmax_) distrmax_=d;
			}
		}
		distrread_=true;
	}
	double statisticalenergy(double phi,double psi) const {
		return -log(distr(phi,psi));
	}
	double statisticalenergy(std::pair<unsigned int,unsigned int> index) const {
			return -log(distr(index)*(double)(phidimension_*psidimension_));
		}
	double phi_statisticalenergy(double phi) const {
		return -log(mdistr_phi(phi));
	}
	double psi_statisticalenergy(double psi) const {
		return -log(mdistr_psi(psi));
	}
	double itplenergy(double phi,double psi,double *dedphi, double *dedpsi) const {
		double hstep=step_*0.5;
		double ecut=0.0;
		std::pair<unsigned int,unsigned int> bottomleft;
		std::pair<unsigned int,unsigned int> topright;
		std::pair<double,double> xy=getcrdingrid(phi,psi,&bottomleft,&topright);
		double e00=statisticalenergy(bottomleft);
		if(e00>ecut) e00=ecut;
		double e10=statisticalenergy(std::make_pair(topright.first,bottomleft.second));
		if(e10>ecut) e10=ecut;
		double e01=statisticalenergy(bottomleft.first,topright.second);
		if(e01>ecut) e01=ecut;
		double e11=statisticalenergy(topright);
		if(e11>ecut) e11=ecut;
		double x=xy.first;
		double y=xy.second;
		assert(x>-0.00000001 && x<1.0000001);
		assert(y>-0.00000001 && y<1.0000001);
		double ene=e00*(1.0-x)*(1-y) + e10*x*(1.0-y) + e01*(1.0-x)*y + e11*x*y;
		*dedphi=((e10-e00)*(1.0-y)+(e11-e01)*y)/step_;
		*dedpsi=((e01-e00)*(1.0-x)+(e11-e10)*x)/step_;
		return ene;
	}
	double intplene_phi(double phi, double *dedphi) const {
		double hstep=0.5*step_;
		double phih=phi-hstep;
		while(phih<0.0) phih+=360.0;
		double ecut=0.0;
		double p0=step_*((double)((int)((phih)/step_))+0.5);
		double e0=phi_statisticalenergy(p0);
		if(e0>ecut) e0=ecut;
		double p1=p0+step_;
		double e1=phi_statisticalenergy(p1);
		if(e1>ecut) e1=ecut;
		double dp1=p1-phi;
		if(dp1<-180.0) dp1+=360.0;
		if(dp1>180.0) dp1-=360.0;
		double x=1.0-dp1/step_;
//		assert(x>-0.0000001 && x<1.00000001);
		if(!(x>-0.0000001 && x<1.00000001)){
//			std::cout <<"phi: " <<phi<<" x:" <<x <<std::endl;
//			exit(1);
			return 0.0;
		}
		double ene=e0*(1.0-x) +e1*x;
		*dedphi=(e1-e0)/step_;
		return ene;
	}
	double intplene_psi(double psi, double *dedpsi) const {
		double hstep=0.5*step_;
		double ecut=0.0;
		double psih=psi-hstep;
		while(psih<0.0) psih+=360.0;
		double p0=step_*((double)((int)((psih)/step_))+0.5);
		double e0=psi_statisticalenergy(p0);
		if(e0>ecut) e0=ecut;
		double p1=p0+step_;
		double e1=psi_statisticalenergy(p1);
		if(e1>ecut) e1=ecut;
		p1=step_*((double)((int)(p1/step_))+0.5);
		double dp1=p1-psi;
		if(dp1>180.0) dp1-=360.0;
		if(dp1<-180.0) dp1+=360.0;
		double x=1.0-dp1/step_;
//		assert(x>-0.0000001 && x<1.00000001);
		if(!(x>-0.0000001 && x<1.00000001)){
//			std::cout <<"phi: " <<phi<<" x:" <<x <<std::endl;
//			exit(1);
			return 0.0;
		}
		double ene=e0*(1.0-x) +e1*x;
		*dedpsi=(e1-e0)/step_;
		return ene;
	}

	template<typename RNG>
	void randomphipsi(RNG &rng,double *phi, double *psi) const {
		while (true) {
			*phi=rng()*360.0;
			*psi=rng()*360.0;
			double dis=distr(*phi,*psi)/(double)(phidimension_*psidimension_);
			if(rng()<dis/distrmax_) return;
		}
	}
	template<typename RNG>
	void randompsi(RNG &rng,double phi, double *psi) const {
		while (true) {
			*psi=rng()*360.0;
			double dis=distr(phi,*psi)/(double)(phidimension_*psidimension_);
			if(rng()<dis/distrmax_) return;
		}
	}
	template<typename RNG>
	void randomphi(RNG &rng,double *phi, double psi) const {
		while (true) {
			*phi=rng()*360.0;
			double dis=distr(*phi,psi)/(double)(phidimension_*psidimension_);
			if(rng()<dis/distrmax_) return;
		}
	}

private:
	double step_;
	unsigned int phidimension_;
	unsigned int psidimension_;
	bool distrread_{false};
	std::vector<double> distr_;
	double distrmax_;
	std::pair<unsigned int,unsigned int> getindex(double phi,double psi) const {
		while(phi<0.0) phi +=360.0;
		while (phi >=360.0) phi -=360.0;
		while (psi<0.0) psi +=360.0;
		while (psi >=360.0) psi -=360.0;
		unsigned int phi_idx=phi/step_;
		unsigned int psi_idx=psi/step_;
		return std::make_pair(phi_idx,psi_idx);
	}
	std::pair<double,double> getphipsi(unsigned int phi_idx, unsigned int psi_idx) const {
		double phi=((double) phi_idx +0.5)*step_;
		double psi=((double) psi_idx+0.5)*step_;
		return std::make_pair(phi,psi);
	}
	std::pair<double,double> getcrdingrid(double phi,double psi,
			std::pair<unsigned int,unsigned int> *bottomleft,std::pair<unsigned int,unsigned int> *topright) const {
		double hstep=0.5*step_;
		*bottomleft=getindex(phi-hstep,psi-hstep);
		*topright=getindex(phi+hstep,psi+hstep);
		std::pair<double,double> orig=getphipsi(bottomleft->first,bottomleft->second);
		double dphi=phi-orig.first;
		double dpsi=psi-orig.second;
		while(dphi>360.0) dphi-=360.0;
		while(dphi<0.0) dphi+=360.0;
		while(dpsi>360.0) dpsi-=360.0;
		while(dpsi<0.0) dpsi+=360.0;
		return std::make_pair(dphi/step_,dpsi/step_);
	}
};


}



#endif /* PDBSTATISTICS_PHIPSIDISTR_H_ */

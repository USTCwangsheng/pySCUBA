/*
 * rminsef.cpp
 *
 *  Created on: 2016年6月8日
 *      Author: hyliu
 */
#include "pdbstatistics/tetrabase.h"
#include<sstream>
using namespace NSPproteinrep;
using namespace NSPgeometry;
using namespace NSPpdbstatistics;
BBD2Matrix::BBD2Matrix(BackBoneSite &s1, BackBoneSite &s2){
		std::vector<XYZ> crd1;
		std::vector<XYZ> crd2;
		double theta=109.5*3.14159265/180.0;
		double t=-120*3.14159265/180.0;
		XYZ n1=s1.getcrd(BackBoneSite::NCRD);
		XYZ c1=s1.getcrd(BackBoneSite::CCRD);
		crd1.push_back(s1.getcrd(BackBoneSite::NCRD));
		crd1.push_back(s1.getcrd(BackBoneSite::CCRD));
		crd1.push_back(s1.cbcrd());
		crd1.push_back(InternaltoXYZ(s1.getcrd(BackBoneSite::CACRD),c1,n1,
				1.5,theta,t));
		crd2.push_back(s2.getcrd(BackBoneSite::NCRD));
		crd2.push_back(s2.getcrd(BackBoneSite::CCRD));
		crd2.push_back(s2.cbcrd());
		crd2.push_back(InternaltoXYZ(s2.getcrd(BackBoneSite::CACRD),crd2[C],crd2[N],
				1.5,theta,t));
		for(auto & x1:crd1) {
			for(auto &x2:crd2) {
				dist2_.push_back(distance2(x1,x2));
			}
		}
		double dmin=10000000000.0;
		double dmax=-1.0;
		for(int i=0; i<dist2_.size();++i) {
			if(dist2_[i] > dmax) {
				dmax=dist2_[i];
				max_index_=i;
			}
			if(dist2_[i] < dmin) {
				dmin=dist2_[i];
				min_index_=i;
			}
		}
		std::pair<int,int> min_pair=indextopair(min_index_);
		std::pair<int,int> max_pair=indextopair(max_index_);
		std::vector<int> atoms1;
		std::vector<int> atoms2;
		for(int i=0;i<4; i++) {
			if(i==min_pair.first || i==max_pair.first) continue;
			atoms1.push_back(i);
		}
		for(int i=0;i<4; i++) {
			if(i==min_pair.second || i==max_pair.second) continue;
			atoms2.push_back(i);
		}
		XYZ center1=0.5*(crd1[atoms1[0]]+crd1[atoms1[1]]);
		XYZ center2=0.5*(crd2[atoms2[0]]+crd2[atoms2[1]]);
		phi_=torsion(crd1[atoms1[0]],center1,center2,crd2[atoms2[0]]);
}
std::pair<int,int> BBD2Matrix::indextopair(int idx) {
	int a=idx/4;
	int b=idx-4*a;
	return std::make_pair(a,b);
}

double EneTables::energy(const Table &tab, double r) {
		const std::vector<double> & crds=tab.r;
		const std::vector<double> & values=tab.ene;
//		if(r>crds.back()) r=crds.back();
		if(r>crds.back()) return 0.0;
		int idx=getindex(crds,r);
		double alpha=(r-crds[idx])/(crds[idx+1]-crds[idx]);
		return values[idx]*alpha +values[idx+1]*(1-alpha);
}

double EneTables::subenergy(const Table &rtab,const std::vector<Table> & ptabs,double r,double p) {
		const std::vector<double> & crds=rtab.r;
		if(r>crds.back()) return 0.0;
		int idx=getindex(crds,r);
		double alpha=(r-crds[idx])/(crds[idx+1]-crds[idx]);
		if(alpha >0.5) idx +=1;
		return energy(ptabs[idx],p);
}

int EneTables::intervale(const std::vector<double> &crds,double c, int lower, int upper) {
		if(upper <= lower+1) return lower;
		int middle=lower +(upper-lower)/2;
		if(c < crds[middle]) upper=middle;
		else lower=middle;
		return intervale(crds,c,lower,upper);
	}
/*
void EneTables::read(const std::string & filename){
	std::ifstream ifs;
	ifs.open(filename.c_str());
	std::string sstype;
	int imin;
	int imax;
	double ratio;
	while(ifs.good()) {
		ifs >> sstype;
		if(ifs.eof()) break;
		std::cout <<"Reading tables of relative distribution of minimum distances for " <<sstype <<std::endl;
		double eref=0.0;
		double rtot=0.0;
		std::vector<EneType> tabletypes;
		for (int m=0; m<144; ++m) {
			ifs >> imin >>imax >>ratio;
			     EneType tabletype=std::make_pair(sstype,std::make_pair(imin,imax));
			     tabletypes.push_back(tabletype);
			     tables_.insert(std::make_pair(tabletype,EneTables::Table()));
			std::vector<double> & rlist=tables_.at(tabletype).r;
			std::vector<double> & vlist=tables_.at(tabletype).ene;
			for (int i=0;i<46;++i) {
				double r;
				double p;
				ifs >>r >>p;
				rlist.push_back(r);
				if (p < 1e-20) p=1e-20;
				p=-log(p);
				vlist.push_back(p);
			} //one orientation type
			rtot +=ratio;
			eref += ratio*vlist.back();
		} //all orientation types
		eref=eref/rtot;
		std::cout <<"Energy shift applied for " << sstype <<" " << eref/rtot <<" " <<rtot <<std::endl;
		for(int m=0;m<144;++m) {
			std::vector<double> & vlist=tables_.at(tabletypes.at(m)).ene;
			for (int i=0; i<46; ++i) {
				vlist.at(i) -= eref;
			}
		}
	} //entype file
}
*/

void EneTables::readwithsubtables(const std::string & filename){
	std::ifstream ifs;
	ifs.open(filename.c_str());
	std::string sstype;
	int imin;
	int imax;
	double ratio;
	int rpoints,phipoints, phieneon;
	char buffer[120];
	ifs.getline(buffer,120);
	std::stringstream ss(buffer);
	std::string s1,s2,s3;
	ss >> s1 >>s2 >>s3;
	rpoints=std::stoi(s1);
	phipoints=std::stoi(s2);
	if(std::stoi(s3) != 0 ) phiene_=true;
	while(ifs.good()) {
		ifs >> sstype;
		if(ifs.eof()) break;
		std::cout <<"Reading tables and subtables of relative distribution of minimum distances for " <<sstype <<std::endl;
		double eref=0.0;
		double rtot=0.0;
		std::vector<EneType> tabletypes;
		for (int m=0; m<144; ++m) {
			ifs >> imin >>imax >>ratio;
			     EneType tabletype=std::make_pair(sstype,std::make_pair(imin,imax));
			     tabletypes.push_back(tabletype);
			     tables_.insert(std::make_pair(tabletype,EneTables::Table()));
			     subtables_.insert(std::make_pair(tabletype,std::vector<Table>(46,Table())));
			std::vector<double> & rlist=tables_.at(tabletype).r;
			std::vector<double> & vlist=tables_.at(tabletype).ene;
			std::vector<Table> & stabs=subtables_.at(tabletype);
			for (int i=0;i<rpoints;++i) {
				double r;
				double p;
				ifs >>r >>p;
				rlist.push_back(r);
				if (p < 1e-20) p=1e-20;
				p=-log(p);
				vlist.push_back(p);
				double phi;
				double pphi;
				std::vector<double> &philist=stabs[i].r;
				std::vector<double> &ephilist=stabs[i].ene;
				for(int j=0; j<phipoints; ++j) {
					ifs >>phi >>pphi;
					philist.push_back(phi);
					ephilist.push_back(-log(pphi));
				}
			} //one orientation type
			rtot +=ratio;
			eref += ratio*vlist.back();
		} //all orientation types
		eref=eref/rtot;
		std::cout <<"Energy shift applied for " << sstype <<" " << eref/rtot <<" " <<rtot <<std::endl;
		for(int m=0;m<144;++m) {
			std::vector<double> & vlist=tables_.at(tabletypes.at(m)).ene;
			for (int i=0; i<rpoints; ++i) {
				vlist.at(i) -= eref;
			}
		}
	} //entype file
}

double TetraBASE::twobody(BackBoneSite &s1, BackBoneSite &s2){
	BBD2Matrix d2m(s1,s2);
	EneTables::EneType et=enetype(s1.sscodechar(),s2.sscodechar(),d2m);
	double rmin=sqrt(d2m.matrix()[d2m.min_index()]);
	double phi=d2m.phi();
	return enetables_.energy(et,rmin,phi);
}

typename EneTables::EneType TetraBASE::enetype(char s1, char s2, const BBD2Matrix & d2m) {
		typename EneTables::OrienType ot=d2m.minmaxindices();
		typename EneTables::SSType ss;
		ss.push_back(s1);
		ss.push_back(s2);
		return std::make_pair(ss,ot);
}


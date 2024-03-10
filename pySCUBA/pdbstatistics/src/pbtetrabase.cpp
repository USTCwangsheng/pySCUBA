/*
 * pbtetrabase.cpp
 *
 *  Created on: 2017年8月4日
 *      Author: hyliu
 */
#include "pdbstatistics/pbtetrabase.h"
#include "dataio/datapaths.h"
#include<sstream>

using namespace NSPpdbstatistics;
using namespace NSPproteinrep;
using namespace NSPgeometry;

double NSPpdbstatistics::tetrapairenergy(const NSPproteinrep::BackBoneSite &s1, const NSPproteinrep::BackBoneSite &s2){
//	static TetraPairETables & etables=gettetrapairetables();
//	static TetraPairETables & etables_hb=gettetrapairetables_hb();
//	static TetraPairETables & etables_hb=gettetrapairetables_nohb();
	TetraGeom geom(s1,s2);
	TetraPairETables &etables=gettetrapairetables(geom.hbond());
	return etables.energy(std::make_pair(geom.pbtypes(),geom.orient()),geom.rmin());
}

/*TetraPairETables & NSPpdbstatistics::gettetrapairetables(const std::string &filename){
	static TetraPairETables etables;
	static bool tableread{false};
	if(!tableread || !filename.empty()){
		std::string fn=filename;
		if(fn.empty()) fn=NSPdataio::datapath()+"pbtetrasef.dat";
		if(!readetables(fn,&etables)){
			std::cout <<"Read tetrapairtable failed. Is the file in place?"<<std::endl;
			abort();
		}
	}
	return etables;
}*/


TetraPairETables & NSPpdbstatistics::gettetrapairetables(bool hbond,const std::string &filename){
	static TetraPairETables etables;
	static TetraPairETables etables_hb;
	static TetraPairETables etables_nohb;

	static bool tableread{false};
	if(!tableread || !filename.empty()){
		std::string fn=filename;
		if(fn.empty()) fn=NSPdataio::datapath()+"pbtetrasef.dat";
		if(!readetables(fn,&etables,&etables_hb,&etables_nohb)){
			std::cout <<"Read tetrapairtable failed. Is the file in place?"<<std::endl;
			abort();
		}
		tableread=true;
	}
	if(hbond) return etables_hb;
	else return etables_nohb;
}
bool NSPpdbstatistics::readetables(const std::string &filename, TetraPairETables *etables,
		TetraPairETables *etables_hb, TetraPairETables *etables_nohb){
	std::ifstream ifs;
	ifs.open(filename.c_str());
	char buffer[120];
	while(ifs.getline(buffer,120)){
		std::string line(buffer);
		std::stringstream sstr(line);
		char pba,pbb;
		int rminidx, rmaxidx;
		double po1,po2;
		sstr>>pba >>pbb >>rminidx >>rmaxidx >>po1 >>po2;
		std::vector<double> r;
		std::vector<double> val;
		std::vector<double> vhb;
		std::vector<double> vnohb;
		for (int i=0;i<45;++i){
			ifs.getline(buffer,120);
			std::string line1(buffer);
			std::stringstream sstr1(line1);
			double ri,vi, hb_cor,nohb_cor;
			sstr1 >>ri >>vi >>hb_cor >>nohb_cor;
			vi *=(po2*144.0);
			hb_cor *=vi;
			nohb_cor *=vi;
			if(vi<1.e-20) vi=1.e-20;
			if(hb_cor <1.e-20) hb_cor=1.e-20;
			if(nohb_cor <1.e-20) nohb_cor=1.e-20;
			r.push_back(ri);
			val.push_back(-log(vi));
			vhb.push_back(-log(hb_cor));
			vnohb.push_back(-log(nohb_cor));
		}
		TetraPairType pt=std::make_pair(std::make_pair(pba,pbb),std::make_pair(rminidx,rmaxidx));
		etables->inserttable(pt,r,val);
		etables_hb->inserttable(pt,r,vhb);
		etables_nohb->inserttable(pt,r,vnohb);
	}
	if(etables->empty()) return false;
	return true;
}
void TetraGeom::dist2matrix(const NSPproteinrep::BackBoneSite &s1,
		const NSPproteinrep::BackBoneSite &s2, std::vector<double> *matrix) {
	{
		matrix->clear();
		std::vector<NSPgeometry::XYZ> crd1, crd2;
		static const double theta = 109.5 * 3.14159265 / 180.0;
		static const double t = -120 * 3.14159265 / 180.0;
		XYZ n1 = s1.getcrd(BackBoneSite::NCRD);
		XYZ c1 = s1.getcrd(BackBoneSite::CCRD);
		crd1.push_back(s1.getcrd(BackBoneSite::NCRD));
		crd1.push_back(s1.getcrd(BackBoneSite::CCRD));
		crd1.push_back(s1.cbcrd());
		crd1.push_back(
				InternaltoXYZ(s1.getcrd(BackBoneSite::CACRD), c1, n1, 1.5,
						theta, t));
		crd2.push_back(s2.getcrd(BackBoneSite::NCRD));
		crd2.push_back(s2.getcrd(BackBoneSite::CCRD));
		crd2.push_back(s2.cbcrd());
		crd2.push_back(
				InternaltoXYZ(s2.getcrd(BackBoneSite::CACRD), crd2[1], crd2[0],
						1.5, theta, t));
		for (auto & r1 : crd1)
			for (auto &r2 : crd2)
				matrix->push_back((r1 - r2).squarednorm());
	}
}
void TetraGeom::init(const std::vector<double> & matrix) {
	int minidx;
	int maxidx;
	double rmin2 = 10000000.0;
	double rmax2 = -1.0;
	int idx = 0;
	for (auto &r2 : matrix) {
		if (rmin2 > r2) {
			minidx = idx;
			rmin2 = r2;
		}
		if (rmax2 < r2) {
			rmax2 = r2;
			maxidx = idx;
		}
		++idx;
	}
	orient_ = std::make_pair(minidx, maxidx);
	rmin_ = sqrt(rmin2);
}

/*
 * estimatepbtetrabase.cpp
 *
 *  Created on: 2017年8月4日
 *      Author: hyliu
 */
#include <backbone/backbonesite.h>
#include <pdbstatistics/pbtetrabase.h>
#include <pdbstatistics/proteinblock.h>
#include <dstl/randomengine.h>
#include <geometry/rotation.h>

#include <iostream>
#include <fstream>
#include <map>
using namespace NSPpdbstatistics;
using namespace NSPproteinrep;
std::map<TetraGeom::PbTypes, std::map<TetraGeom::OrientType, double>> orientcount;
std::map<TetraGeom::PbTypes, std::map<TetraGeom::OrientType, double>> shell_orientcount;
std::map<TetraGeom::PbTypes,double> pbtypescount_total;
std::map<TetraGeom::PbTypes,double> shell_pbtypescount;
std::map<TetraGeom::OrientType,std::vector<double>> bghists;
std::map<TetraGeom::OrientType,std::vector<double>> hbfraction_bg;
std::map<TetraGeom::PbTypes,std::map<TetraGeom::OrientType,std::vector<double>>> rminhist;
std::map<TetraGeom::PbTypes,std::map<TetraGeom::OrientType,std::vector<double>>> hbfraction;
std::vector<double> refhist;
std::map<TetraGeom::PbTypes,std::map<TetraGeom::OrientType,std::vector<double>>> refhbfraction;
std::map<TetraGeom::PbTypes,std::map<TetraGeom::OrientType,std::vector<double>>> refhbhist;
std::vector<HBondGeometry> hbgs;
double rstep { 0.25 };
double npseudo{100.0};
long paircount{1};
void addhistgram(TetraGeom & geom) {
	double rmin = geom.rmin();
	++paircount;
	auto pbs = geom.pbtypes();
	auto ort = geom.orient();
	int ridx = rmin / rstep;
	if (ridx > 44)
		ridx = 44;
	if (orientcount.find(pbs) == orientcount.end())
		orientcount.insert(
				std::make_pair(pbs, std::map<TetraGeom::OrientType, double>()));
	if (orientcount.at(pbs).find(ort) == orientcount.at(pbs).end())
		orientcount.at(pbs).insert(std::make_pair(ort, 0.0));
	if (rminhist.find(pbs) == rminhist.end())
		rminhist.insert(std::make_pair(pbs, std::map<TetraGeom::OrientType,std::vector<double>>()));
	if(rminhist[pbs].find(ort) == rminhist[pbs].end())
		rminhist[pbs].insert(std::make_pair(ort,std::vector<double>(45,0.0)));
	if(bghists.find(ort) == bghists.end()) {
		bghists.insert(std::make_pair(ort,std::vector<double>(45,0.0)));
	}
	orientcount[pbs][ort] += 1.0;
	rminhist[pbs][ort][ridx] += 1.0;
	if(!(pbs.first == 'm' || pbs.first =='d' || pbs.second =='m' ||pbs.second =='d')){
		bghists[ort][ridx] +=1.0;
	}
	if(shell_orientcount.find(pbs) == shell_orientcount.end()){
		shell_orientcount.insert(std::make_pair(pbs,std::map<TetraGeom::OrientType,double>()));
	}
	if(shell_orientcount[pbs].find(ort)==shell_orientcount[pbs].end()){
		shell_orientcount[pbs].insert(std::make_pair(ort,0.0));
	}
	if(shell_pbtypescount.find(pbs) ==shell_pbtypescount.end())
		shell_pbtypescount.insert(std::make_pair(pbs,0.0));
	if(ridx==44) {
		shell_orientcount[pbs][ort] +=1.0;
		shell_pbtypescount[pbs] +=1.0;
	}
	if(hbfraction_bg.find(ort) == hbfraction_bg.end()){
		hbfraction_bg.insert(std::make_pair(ort,std::vector<double>(45,0.0)));
	}
	if(hbfraction.find(pbs) == hbfraction.end()) {
		hbfraction.insert(std::make_pair(pbs,std::map<TetraGeom::OrientType,std::vector<double>>()));
	}
	if(hbfraction[pbs].find(ort)==hbfraction[pbs].end()){
		hbfraction[pbs].insert(std::make_pair(ort,std::vector<double>(45,0.0)));
	}
	if(geom.hbond()) {
		if(!(pbs.first == 'm' || pbs.first =='d' || pbs.second =='m' ||pbs.second =='d')){
			hbfraction_bg[ort][ridx] +=1.0;
		}
		hbfraction[pbs][ort][ridx] +=1.0;
		if(hbgs.size()<20000) hbgs.push_back(geom.hbbondgeometry());
	}
}
void normalize() {
	for(auto & ordis:bghists){
		for(int i=0;i<45;++i){
			double rdis=ordis.second[i];
			if(rdis>0) hbfraction_bg[ordis.first][i] /=rdis;
			else hbfraction_bg[ordis.first][i]=0.0;
		}
		double totalcount=0.0;
		for(auto &rdis:ordis.second){
			totalcount +=rdis;
		}
		int ridx=0;
		for(auto &rdis:ordis.second) {
			rdis= npseudo*(rdis/totalcount);
			hbfraction_bg[ordis.first][ridx++] *=rdis;
		}
		std::cout <<ordis.first.first << " "<<ordis.first.second<<" "
				<<ordis.second.back()<<std::endl;
	}
	for(auto &ortdistr:orientcount){
		for(auto &ort:ortdistr.second)
			ort.second +=0.05*npseudo;
		for(auto &ort:shell_orientcount[ortdistr.first]){
			ort.second +=0.05*bghists[ort.first].back();
			std::cout <<ortdistr.first.first << " "<<ortdistr.first.second<<" "
					<<ort.first.first<<" "<<ort.first.second<<" "
					<<ort.second<<std::endl;
		}
	}
	for (auto &ortdistr : orientcount) {
		double totalcount = 0.0;
		for (auto &ort : ortdistr.second)
			totalcount += ort.second;
		double shell_count=0.0;
		for(auto  &ort:shell_orientcount[ortdistr.first])
			shell_count +=ort.second;
		std::cout << ortdistr.first.first <<" "<<ortdistr.first.second <<" "
				<< shell_count/
				(pbtypescount_total[ortdistr.first]-totalcount)<<" " <<shell_count <<" "
				<< totalcount <<" "<< pbtypescount_total[ortdistr.first]<<std::endl;
		for (auto &ort : ortdistr.second)
			ort.second /= totalcount;
		for(auto &ort:shell_orientcount[ortdistr.first])
			ort.second /=shell_count;
	}

	for (auto &rmindistr : rminhist) {
		for (auto &ordis : rmindistr.second) {
			auto pbtypes=rmindistr.first;
			auto orient=ordis.first;
			double totalcount=0.0;
			int ridx=0;
			for(auto &rdis:ordis.second){
				totalcount += rdis;
				rdis +=bghists[ordis.first][ridx++];
			}
			for(int i=0;i<45;++i) {
				hbfraction[pbtypes][orient][i]+= hbfraction_bg[orient][i];
				if(ordis.second[i] >0)
					hbfraction[pbtypes][orient][i] /=ordis.second[i];
				else
					hbfraction[pbtypes][orient][i] =0.0;
			}
			std::cout <<pbtypes.first<<" " <<pbtypes.second <<" " << orient.first <<" " <<orient.second;
			std::cout <<" " <<totalcount<<std::endl;
			totalcount += npseudo;
			for (auto &rdis : ordis.second)
				rdis /= totalcount;
		}
	}

}

void writeresult(std::ostream &os) {
	for(auto &rmindist:rminhist){
		auto pbtypes=rmindist.first;
		for(auto &ordis:rmindist.second){
			auto orient=ordis.first;
			os <<pbtypes.first<<" " <<pbtypes.second <<" " << orient.first <<" " <<orient.second;
			os <<" "<<orientcount[pbtypes][orient]<<" "<<shell_orientcount[pbtypes][orient]<<std::endl;
			for(int i=0;i<45;++i) {
				ordis.second[i] /=refhist[i];
			}
			double nc=ordis.second.back();
			for(int i=0;i<45;++i) {
				double r=rstep*(double) i + 0.5*rstep;
				double res1=1.0;
				double res2=1.0;
				if(hbfraction[pbtypes][orient][i] >0) {
				   double deno=refhbfraction[pbtypes][orient][i];
				   if(deno==0) deno=1.0/(double) refhbhist[pbtypes][orient][i];
				   res1=hbfraction[pbtypes][orient][i]/deno;
				   res2=(1-hbfraction[pbtypes][orient][i])/(1-deno);
				}
				os << r <<" "<<ordis.second[i]/nc<<" ";
				os <<res1<<" "<<res2<<std::endl;
			}
		}
	}
}
void makerefhist(const std::vector<BackBoneSite> &sites) {
	auto & reng = NSPdstl::RandomEngine<>::getinstance();
	refhist.resize(45, 0.0);
	for (int i = 0; i < 10000; ++i) {
		int posi1;
		char pbtype1;
		while (true) {
			posi1=reng.intrng(2, sites.size() - 3)();
			if (!fragstartsite(sites.begin()+posi1 - 2, sites.end(), 5)) continue;
			pbtype1=ProteinBlock::pbtype(sites, posi1);
			if(pbtype1 =='m' ||pbtype1=='d') continue;
			break;
		}
		int posi2;
		char pbtype2;
		while (true) {
			posi2=reng.intrng(2, sites.size() - 3)();
			if (!fragstartsite(sites.begin()+posi2 - 2, sites.end(), 5)) continue;
			pbtype2=ProteinBlock::pbtype(sites, posi2);
			if(pbtype2 =='m' ||pbtype2=='d') continue;
			break;
		}
		BackBoneSite s1 = sites.at(posi1);
		BackBoneSite s2=  sites.at(posi2);
		s1.sscode = pbtype1;
		s2.sscode = pbtype2;
		s1.translate(-1.0 * s1.cacrd());
		for (int n = 0; n < 50; ++n) {
			s2.translate(-1.0* s2.cacrd());
			NSPgeometry::XYZ trans(reng.realrng(0.0, 1.0), 1.0);
			trans = 14.0 * trans;
			s2.translate(trans);
			NSPgeometry::XYZ axis(reng.realrng(), 1.0);
			double angle = reng.realrng()() * 180.0;
			NSPgeometry::Rotation rot;
			rot.init(NSPgeometry::QuaternionCrd(axis, angle),
					s2.cacrd());
			s2.rotate(rot);
			TetraGeom geom(s1, s2);
			double r = geom.rmin();
			int ridx = r / rstep;
			if (ridx > 44)
				ridx = 44;
			refhist[ridx] += 1.0;
		}
	}
	double totalcount = 0.0;
	for (auto &c : refhist)
		totalcount += c;
	for (auto &c : refhist)
		c /= totalcount;
}
void makehbref(const std::vector<BackBoneSite> &sites) {
	auto & reng = NSPdstl::RandomEngine<>::getinstance();
	for (int i = 0; i < 100000; ++i) {
		int posi1;
		char pbtype1;
		while (true) {
			posi1=reng.intrng(2, sites.size() - 3)();
			if (!fragstartsite(sites.begin()+posi1 - 2, sites.end(), 5)) continue;
			pbtype1=ProteinBlock::pbtype(sites, posi1);
			if(pbtype1 =='m' ||pbtype1=='d') {
				if(reng.realrng(0.0,1.0)()>0.2)continue;
			}
			break;
		}
		int posi2;
		char pbtype2;
		while (true) {
			posi2=reng.intrng(2, sites.size() - 3)();
			if (!fragstartsite(sites.begin()+posi2 - 2, sites.end(), 5)) continue;
			pbtype2=ProteinBlock::pbtype(sites, posi2);
			if(pbtype2 =='m' ||pbtype2=='d') {
					if(reng.realrng(0.0,1.0)()>0.2)continue;
				}
			break;
		}
		BackBoneSite s1 = sites.at(posi1);
		BackBoneSite s2=  sites.at(posi2);
		s1.sscode = pbtype1;
		s2.sscode = pbtype2;
		int ntry=5000;
		if(pbtype1=='j' || pbtype1=='e' || pbtype1=='n' ||pbtype1=='i') ntry=20000;
		if(pbtype2=='j' || pbtype2=='e' || pbtype2=='n' ||pbtype2=='i') ntry=20000;
		s1.translate(-1.0 * s1.cacrd());
		for (int n = 0; n < ntry; ++n) {
			s2.translate(-1.0* s2.cacrd());
			NSPgeometry::XYZ trans(reng.realrng(0.0, 1.0), 1.0);
			double norm=sqrt(trans.squarednorm());
			trans = (reng.realrng(0.375,1.0)()*8.0/norm) * trans;
			s2.translate(trans);
			NSPgeometry::XYZ axis(reng.realrng(), 1.0);
			double angle = reng.realrng(0.0,1.0)() * 180.0;
			NSPgeometry::Rotation rot;
			rot.init(NSPgeometry::QuaternionCrd(axis, angle),
					s2.cacrd());
			s2.rotate(rot);
			TetraGeom geom(s1, s2);
			double r = geom.rmin();
			int ridx = r / rstep;
			if (ridx > 44)
				ridx = 44;
			auto ort=geom.orient();
			auto pbs=geom.pbtypes();
			if(refhbhist.find(pbs)==refhbhist.end()){
				refhbhist.insert(std::make_pair(pbs,std::map<TetraGeom::OrientType,std::vector<double>>()));
				refhbfraction.insert(std::make_pair(pbs,std::map<TetraGeom::OrientType,std::vector<double>>()));
			}
			if(refhbhist[pbs].find(ort)==refhbhist[pbs].end()){
				refhbhist[pbs].insert(std::make_pair(ort,std::vector<double>(45,0.0)));
				refhbfraction[pbs].insert(std::make_pair(ort,std::vector<double>(45,0.0)));
			}
			refhbhist[pbs][ort][ridx] +=1.0;
			if(geom.hbond())refhbfraction[pbs][ort][ridx] +=1.0;
		}
	}
	for(auto &refohis:refhbhist){
		for(auto &ohis:refohis.second){
			for(int i=0;i<45;++i) {
				std::cout <<"hbref "<<refohis.first.first <<" "<< refohis.first.second
						<<" " <<ohis.first.first<<" "<<ohis.first.second<<" "
						<< ohis.second[i] <<" "<< refhbfraction[refohis.first][ohis.first][i] <<std::endl;
				if(ohis.second[i] >0) refhbfraction[refohis.first][ohis.first][i] /=ohis.second[i];
				else refhbfraction[refohis.first][ohis.first][i]=0.0;
			}
		}
	}
}
int main(int argc, char **argv) {
	std::vector<BackBoneSite> sites;
	readbackbonesites(std::string(argv[1]), sites);
	int minsep = 6;
	makerefhist(sites);
	makehbref(sites);
	auto start2 = sites.begin() + 2;
	int nchain=0;
	for (auto iter1 = start2 + minsep; iter1 != sites.end()-2; ++iter1) {
		if (chainstartsite(iter1)) {
			std::cout<<"Number of chains: " <<++nchain <<" Number of Pairs: " <<paircount <<std::endl;
			start2=iter1+2;
		}
		if(iter1-start2 < minsep) continue;
		if (!fragstartsite(iter1 - 2, sites.end(), 5))
			continue;
		BackBoneSite s1 = *iter1;
		s1.sscode = ProteinBlock::pbtype(sites, iter1 - sites.begin());
		for (auto iter2 = start2; iter2 != iter1 - minsep + 1; ++iter2) {
			if (!fragstartsite(iter2 - 2, sites.end(), 5))
				continue;
			BackBoneSite s2 = *iter2;
			s2.sscode = ProteinBlock::pbtype(sites, iter2 - sites.begin());
			double r12=(s1.cacrd()-s2.cacrd()).squarednorm();
			TetraGeom geom1(s1, s2);
			if(geom1.rmin()<1.0) {
				std::cout <<s1.toString();
				std::cout<<s2.toString();
			}
			if(pbtypescount_total.find(geom1.pbtypes()) == pbtypescount_total.end())
				pbtypescount_total.insert(std::make_pair(geom1.pbtypes(),0.0));
			pbtypescount_total[geom1.pbtypes()] +=2.0;
			if(r12>196.0) continue;
			addhistgram(geom1);
			TetraGeom geom2(s2, s1);
			addhistgram(geom2);
		}
//		if (hbgs.size()>=20000) break;
	}
	normalize();
	std::ofstream ofs;
	ofs.open("pbtetrasef_hb.dat");
	writeresult(ofs);
	ofs.close();
	ofs.open("hbgeometry.dat");
	for(auto & hbg:hbgs){
		ofs <<hbg.rda <<" "<<hbg.rha <<" " <<hbg.adab <<" "<<hbg.adha <<" "<<hbg.ahab<<std::endl;
	}
}


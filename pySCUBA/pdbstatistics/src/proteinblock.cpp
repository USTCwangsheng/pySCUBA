/*
 * proteinblock.cpp
 *
 *  Created on: 2017年6月13日
 *      Author: hyliu
 */
#include "pdbstatistics/proteinblock.h"
#include <cmath>
using namespace NSPpdbstatistics;
using namespace NSPproteinrep;
std::map<char,std::vector<double>> ProteinBlock::stdpbs
{
	{'a', { 41.14,   75.53,  13.92,  -99.80,  131.88,  -96.27, 122.08,  -99.68}},
	{'b', {108.24,  -90.12, 119.54,  -92.21,  -18.06, -128.93, 147.04,  -99.90}},
	{'c', {-11.61, -105.66,  94.81, -106.09,  133.56, -106.93, 135.97, -100.63}},
	{'d', {141.98, -112.79, 132.20, -114.79,  140.11, -111.05, 139.54, -103.16}},
	{'e', {133.25, -112.37, 137.64, -108.13,  133.00,  -87.30, 120.54,   77.40}},
	{'f', {116.40, -105.53, 129.32,  -96.68,  140.72,  -74.19, -26.65,  -94.51}},
	{'g', {0.40,  -81.83,   4.91, -100.59,   85.50,  -71.65, 130.78,   84.98}},
	{'h', {119.14, -102.58, 130.83,  -67.91,  121.55,   76.25,  -2.95,  -90.88}},
	{'i', {130.68,  -56.92, 119.26,   77.85,   10.42,  -99.43, 141.40,  -98.01}},
	{'j', {114.32, -121.47, 118.14,   82.88, -150.05,  -83.81,  23.35,  -85.82}},
	{'k', {117.16,  -95.41, 140.40,  -59.35,  -29.23,  -72.39, -25.08,  -76.16}},
	{'l', {139.20,  -55.96, -32.70,  -68.51,  -26.09,  -74.44, -22.60,  -71.74}},
	{'m', {-39.62,  -64.73, -39.52,  -65.54,  -38.88,  -66.89, -37.76,  -70.19}},
	{'n', {-35.34,  -65.03, -38.12,  -66.34,  -29.51,  -89.10,  -2.91,   77.90}},
	{'o', {-45.29,  -67.44, -27.72,  -87.27,    5.13,   77.49,  30.71,  -93.23}},
	{'p', {-27.09,  -86.14,   0.30,   59.85,   21.51,  -96.30, 132.67,  -92.91}}
};
char ProteinBlock::pbtype(const std::vector<double> &dihs,double *dev){
	char res=' ';
	double devmin=100000000.0;
	for(auto & pb:stdpbs) {
		double sum2=0.0;
		for(int i=0; i<8;++i) {
			double diff= dihdiff(dihs[i],pb.second[i]);
			sum2 += diff*diff;
		}
		if(sum2 < devmin) {
			devmin=sum2;
			res=pb.first;
		}
	}
	if(dev) *dev=sqrt(devmin/8.0);
	return res;
}

char ProteinBlock::pbtype(const std::vector<NSPproteinrep::BackBoneSite> &sites,
		int posi,double *dev){
	if(posi<2 || posi >=sites.size()-2) return 't';
	std::vector<double> dihs=extractdihvec(sites,posi);
	if(dihs.empty()) return 't';
	return pbtype(dihs,dev);
}
char ProteinBlock::pbtype(const std::vector<NSPproteinrep::BackBoneSite> &sites, bool keepss,
		int posi, double *dev){
	if(keepss){
		char sscode=sites[posi].sscode;
		if(sscode =='H' || sscode=='E')  return sscode;
	}
	return pbtype(sites,posi,dev);
}
std::vector<double> ProteinBlock::deviations(char pb,const std::vector<double> &dihs){
	if(stdpbs.find(pb) == stdpbs.end()) return std::vector<double>();
	std::vector<double> &dihref=stdpbs.at(pb);
	std::vector<double> dev;
	for(int i=0; i<8;++i) {
			dev.push_back(dihdiff(dihs[i],dihref[i]));
	}
	return dev;
}
std::vector<double> ProteinBlock::deviations(char pb, const std::vector<NSPproteinrep::BackBoneSite> &sites,
		int posi){
	std::vector<double> dihs=extractdihvec(sites,posi);
	return deviations(pb,dihs);
}
std::vector<double> ProteinBlock::extractdihvec(const std::vector<NSPproteinrep::BackBoneSite> &sites,
		int posi){
	assert(posi >=2 && posi <sites.size()-2);
	std::vector<double> dihs;
	for(int i=posi-2; i<=posi+2;++i){
		if(sites[i].isgap) return dihs;
	}
	dihs.push_back(sites[posi-2].psi());
	dihs.push_back(sites[posi-1].phi());
	dihs.push_back(sites[posi-1].psi());
	dihs.push_back(sites[posi].phi());
	dihs.push_back(sites[posi].psi());
	dihs.push_back(sites[posi+1].phi());
	dihs.push_back(sites[posi+1].psi());
	dihs.push_back(sites[posi+2].phi());
	return dihs;
}
void ProteinBlock::setpbtypes(std::vector<BackBoneSite> & chain,int posistart,int posiend)
{
		if (posistart<0) posistart=0;
		if(posiend<0||posiend>chain.size()) posiend=chain.size();
		for (int posi=posistart;posi<posiend;++posi){
			chain[posi].sscode=pbtype(chain,posi);
		}
}
void ProteinBlock::setpbtypes(std::vector<BackBoneSite> & chain,bool keepss,int posistart,int posiend)
{
		if (posistart<0) posistart=0;
		if(posiend<0||posiend>chain.size()) posiend=chain.size();
		if(!keepss) {
			setpbtypes(chain,posistart,posiend);
			return;
		} else {
			for (int posi=posistart;posi<posiend;++posi){
				char & ss=chain[posi].sscode;
				if(ss != 'H' && ss != 'E') ss=pbtype(chain,posi);
			}
		}
}

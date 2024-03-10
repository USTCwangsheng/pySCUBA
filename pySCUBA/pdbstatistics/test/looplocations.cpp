/*
 * looplocations.cpp
 *
 *  Created on: 2016年12月10日
 *      Author: hyliu
 */

#include "backbone/backbonesite.h"

#include <iostream>
#include <fstream>
#include <map>


using namespace NSPproteinrep;
int minlooplength = 2;
int maxlooplength = 15;
void loopinchains(std::vector<BackBoneSite> *sites, std::ostream &os) {
	std::vector<long> elembegins; //the starting location of SS elements
	elembegins.push_back(0);
	char ssold = sites->at(0).sscodechar();
	long chainstart = 0;
	long chainend = 0;
	std::map<long, long> chainrange;
	for (auto it = sites->begin() + 1; it != sites->end(); ++it) {
		if (chainendsite(it, sites->end())) {
			chainend = it - sites->begin();
			chainrange.insert(std::make_pair(chainstart, chainend));
		}
		if (chainstartsite(it)) {
			chainstart = it - sites->begin();
			elembegins.push_back(it - sites->begin());
			ssold = it->sscodechar();
		} else {
			if (it->sscodechar() == ssold)
				continue;
			elembegins.push_back(it - sites->begin());
			ssold = it->sscodechar();
		}
	}
	std::map<int,std::pair<double,double>> distranges;
	for(int i=minlooplength;i <=maxlooplength;++i){
		distranges.insert(std::make_pair(i,std::make_pair<double,double>(100000.0,-1.0)));
	}
	elembegins.push_back(sites->size());
	for (auto it = elembegins.begin(); it != elembegins.end() - 1; ++it) {
		if (chainrange.find(*it) != chainrange.end()) {
			chainstart = *it;
			chainend = chainrange[chainstart];
		}
		int length = *(it + 1) - *it;
		if (length < minlooplength || length > maxlooplength)
			continue;
		if (!fragstartsite(sites->begin() + (*it) - 1, sites->end(), length + 2,
				std::string(), false))
			continue;
		BackBoneSite & s = sites->at(*it);
		if (s.sscodechar() != 'C')
			continue;
		else {
			bool start = *it == 0;
			if (!start)
				start = chainstartsite(sites->begin() + (*it));
			if (start)
				continue;
			long pos0 = *(it - 1);
			BackBoneSite &s0 = sites->at(pos0);
			long pos1 = *(it + 1);
			if (chainendsite(sites->begin() + pos1 - 1, sites->end()))
				continue;
			BackBoneSite &n_end = sites->at(*it);
			BackBoneSite &c_end = sites->at(pos1 - 1);
			double dist = NSPgeometry::distance(
					n_end.getcrd(BackBoneSite::CACRD),
					c_end.getcrd(BackBoneSite::CACRD));
			if(dist< distranges[length].first)  distranges[length].first=dist;
			if(dist >distranges[length].second) distranges[length].second=dist;
			os << chainstart << "\t" << *it << "\t" << pos1 << "\t"
					<< chainend + 1 << "\t" << length << "\t" << dist
					<< std::endl;
		}
	} //elementbegins
	std::cout <<"Looplength   MINDIST   MAXDIST" <<std::endl;
	for(int i=minlooplength;i <=maxlooplength;++i){
		std::cout <<"\t"<<i <<"\t" <<distranges[i].first<<"\t"
					<<distranges[i].second<<std::endl;
	}
}

int main(int argc, char **argv) {
	if (argc < 2) {
		std::cout << "usage: " << argv[0] << " sitefile\n";
		exit(0);
	}
	std::vector<BackBoneSite> sites;
	readbackbonesites(std::string(argv[1]), sites);
	std::ofstream ofs;
	ofs.open("loop_locations_lengths_distances.dat");
	loopinchains(&sites, ofs);
}

/*
 * testsselements.cpp
 *
 *  Created on: 2016年5月17日
 *      Author: hyliu
 */
#include "pdbstatistics/sselements.h"
#include "pdbstatistics/tetrabase.h"
#include "dataio/splitstring.h"
using namespace NSPproteinrep;
using namespace NSPpdbstatistics;
void readsegments(std::istream &ifs, std::vector<std::vector<BackBoneSite>> & segments){
	char buffer[120];
	ifs.getline(buffer,120);
	std::vector<int> segmentlengths=NSPdataio::integersInString(std::string(buffer));
	segments.resize(segmentlengths.size(),std::vector<BackBoneSite>());
	for(int i=0;i<segmentlengths.size();++i) {
		if(!readbackbonesites(ifs,segmentlengths[i],segments.at(i))){
			std::cout <<"Error reading segments."<<std::endl;
			exit(1);
		}
	}
}
int main(int argc, char **argv){
	std::vector<BackBoneSite> sites;
	readbackbonesites(std::string(argv[1]), sites);
	SSElements sse(&sites);
	std::vector<std::vector<BackBoneSite>> segments;
	std::ifstream ifs;
	ifs.open(argv[2]);
	readsegments(ifs,segments);
	for(unsigned int i=0; i<segments.size();++i) {
		for(unsigned int j=0; j<segments.size();++j) {
			if(i==j) continue;
			BackBoneSite & head=segments[i].back();
			BackBoneSite & tail=segments[j][0];
			std::map<unsigned int,std::pair<unsigned int,unsigned int>> nmatches;
			BBD2Matrix d2m(head,tail);
			typename EneTables::EneType ltype=TetraBASE::enetype(head.sscodechar(),tail.sscodechar(),d2m);
			std::string looptype(2,' ');
			looptype[0]=head.sscodechar();
			looptype[1]=tail.sscodechar();
			double r2=sqrt(d2m.matrix()[d2m.min_index()]);
			std::vector<std::vector<long>> *templates=sse.loopsptr()[looptype];
			for( auto it = templates->begin(); it != templates->end(); ++it) {
				BackBoneSite &htemp= sites[(*it)[1]-1];
				BackBoneSite &ttemp= sites[(*it)[2]];
				BBD2Matrix td2m(htemp,ttemp);
				typename EneTables::EneType tltype=TetraBASE::enetype(htemp.sscodechar(),ttemp.sscodechar(),td2m);
				if(tltype != ltype) continue;
				double tr2=sqrt(td2m.matrix()[td2m.min_index()]);
				double diff=(r2-tr2)/r2;
				if(diff <-0.1 || diff > 0.1) continue;
				unsigned int looplength=(*it)[2]-(*it)[1];
				if( nmatches.find (looplength) == nmatches.end())
					nmatches.insert(std::make_pair(looplength, std::make_pair(0u,0u)));
				if(diff < 0) nmatches[looplength].first +=1;
				else nmatches[looplength].second +=1;
			}
			std::cout <<"segments: " <<i <<" " <<j <<std::endl;
			for(auto & l:nmatches) {
				std::cout <<"linkerlength " << l.first
						<< " numbers matched: " << l.second.first <<" " <<l.second.second<<std::endl;
			}
		}
	}
	/*
	for (int i=1;i<40; i++) {
		std::cout <<"Length " <<i <<" helices: " <<sse.helices(i).size() <<std::endl;
		for (auto &s:sse.helices(i)){
			BackBoneSite &bs=sites[s];
			std::cout << bs.pdbid <<bs.resid <<std::endl;
		}
		std::cout <<"Length " <<i <<" strands: " <<sse.strands(i).size()<<std::endl;
		for (auto &s:sse.strands(i)){
			BackBoneSite &bs=sites[s];
			std::cout << bs.pdbid <<bs.resid <<std::endl;
		}
	}
	std::cout <<"hhloops: " <<sse.hhloops().size()<<std::endl;
	for(auto &s:sse.hhloops()){
			BackBoneSite &bs=sites[s[1]];
			std::cout <<bs.pdbid<<bs.resid<<"\t"<< s[1]-s[0] <<":"<<s[2]-s[1]<<":"<<s[3]-s[2]<<std::endl;
		}
	std::cout <<"eeloops: "<<sse.eeloops().size() <<std::endl;
	for(auto &s:sse.eeloops()){
		BackBoneSite &bs=sites[s[1]];
		std::cout <<bs.pdbid<<bs.resid<< "\t"<< s[1]-s[0] <<":"<<s[2]-s[1]<<":"<<s[3]-s[2]<<std::endl;
	}
	std::cout <<"ehloops: " <<sse.ehloops().size()<<std::endl;
	for(auto &s:sse.ehloops()){
			BackBoneSite &bs=sites[s[1]];
			std::cout <<bs.pdbid<<bs.resid<< "\t"<<s[1]-s[0] <<":"<<s[2]-s[1]<<":"<<s[3]-s[2]<<std::endl;
		}
	std::cout <<"heloops: " <<sse.heloops().size()<<std::endl;
	for(auto &s:sse.heloops()){
		BackBoneSite &bs=sites[s[1]];
		std::cout <<bs.pdbid<<bs.resid<<"\t"<< s[1]-s[0] <<":"<<s[2]-s[1]<<":"<<s[3]-s[2]<<std::endl;
	}
	std::cout <<"neloops: " <<sse.neloops().size()<<std::endl;
	for(auto &s:sse.neloops()){
		BackBoneSite &bs=sites[s[0]];
		std::cout <<bs.pdbid<<bs.resid<<"\t"<< s[1]-s[0] <<":"<<s[2]-s[1]<<std::endl;
	}
	std::cout <<"nhloops: " <<sse.nhloops().size()<<std::endl;
	for(auto &s:sse.nhloops()){
		BackBoneSite &bs=sites[s[0]];
		std::cout <<bs.pdbid<<bs.resid<< "\t"<<s[1]-s[0] <<":"<<s[2]-s[1]<<std::endl;
	}
	std::cout <<"ecloops: " <<sse.ecloops().size()<<std::endl;
	for(auto &s:sse.ecloops()){
		BackBoneSite &bs=sites[s[1]];
		std::cout <<bs.pdbid<<bs.resid<<"\t"<< s[1]-s[0] <<":"<<s[2]-s[1]<<std::endl;
	}
	std::cout <<"hcloops: " <<sse.hcloops().size()<<std::endl;
	for(auto &s:sse.hcloops()){
		BackBoneSite &bs=sites[s[1]];
		std::cout <<bs.pdbid<<bs.resid<< "\t"<<s[1]-s[0] <<":"<<s[2]-s[1]<<std::endl;
	}
	*/
}


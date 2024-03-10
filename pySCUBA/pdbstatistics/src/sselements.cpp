/*
 * sselements.cpp
 *
 *  Created on: 2016年5月17日
 *      Author: hyliu
 */


#include <pdbstatistics/sselements.h>

using namespace NSPpdbstatistics;
using namespace NSPproteinrep;

SSElements::SSElements(std::vector<BackBoneSite> *sites) {
	sites_=sites;
	std::vector<long> elembegins;
	elembegins.push_back(0);
	char ssold=sites->at(0).sscodechar();
	for (auto it=sites->begin()+1; it != sites->end(); ++it) {
		if(chainstartsite(it)) {
			elembegins.push_back(it-sites->begin());
			ssold=it->sscodechar();
		} else {
			if(it->sscodechar() == ssold) continue;
			elembegins.push_back(it-sites->begin());
			ssold=it->sscodechar();
		}
	}
	elembegins.push_back(sites->size());
	for (auto it=elembegins.begin(); it != elembegins.end()-1; ++it) {
		int length=*(it+1)-*it;
		if(!fragstartsite(sites->begin()+(*it),sites->end(),length,std::string(),false)) continue;
		BackBoneSite & s=sites->at(*it);
		if(s.sscodechar()=='H'){
			if(helices_.find(length) == helices_.end()) helices_.insert(std::make_pair(length,std::vector<long>()));
					helices_[length].push_back(*it);
		}else if(s.sscodechar()=='E') {
			if(strands_.find(length) == strands_.end()) strands_.insert(std::make_pair(length,std::vector<long>()));
			strands_[length].push_back(*it);
		} else {
			bool start= *it==0;
			if( !start) start=chainstartsite(sites->begin()+(*it));
			if(start) {
				long pos1=*(it+1) ;
				if(!chainendsite(sites->begin()+pos1-1,sites->end())) {
					BackBoneSite &s2= sites->at(pos1);
					long pos2=*(it+2);
					if(!fragstartsite(sites->begin()+*it,sites->end(),pos2-*it,std::string(),false)) continue;
					if(s2.sscodechar()=='H') nhloops_.push_back(std::vector<long>({*it,pos1,pos2}));
					else if(s2.sscodechar()=='E') neloops_.push_back(std::vector<long>({*it,pos1,pos2}));
				}
			} else {
				long pos0=*(it-1);
				BackBoneSite &s0=sites->at(pos0);
				long pos1=*(it+1);
				if(!chainendsite(sites->begin()+pos1-1,sites->end())){
					BackBoneSite &s2=sites->at(pos1);
					long pos2=*(it+2);
					if(!fragstartsite(sites->begin()+pos0,sites->end(),pos2-pos0,std::string(),false)) continue;
					if(s0.sscodechar() =='H' && s2.sscodechar() =='H') {
						hhloops_.push_back(std::vector<long>{pos0,*it,pos1,pos2});
					} else if (s0.sscodechar() =='E' && s2.sscodechar() =='E') {
						eeloops_.push_back(std::vector<long>{pos0,*it,pos1,pos2});
					} else if(s0.sscodechar() =='H' && s2.sscodechar() =='E') {
						heloops_.push_back(std::vector<long>{pos0,*it,pos1,pos2});
					}else if (s0.sscodechar() =='E' && s2.sscodechar() =='H') {
						ehloops_.push_back(std::vector<long>{pos0,*it,pos1,pos2});
					}
				} else {
					if(!fragstartsite(sites->begin()+pos0,sites->end(),pos1-pos0,std::string(),false)) continue;
					if(s0.sscodechar()=='H') hcloops_.push_back(std::vector<long>({pos0,*it,pos1}));
					else if(s0.sscodechar()=='E') ecloops_.push_back(std::vector<long>({pos0,*it,pos1}));
				} //chainend
			} // chainstart
		} //sscodechar
	} //elementbegins
	loopsptr_.insert(std::make_pair("HH",&hhloops_));
	loopsptr_.insert(std::make_pair("HE",&heloops_));
	loopsptr_.insert(std::make_pair("EH",&ehloops_));
	loopsptr_.insert(std::make_pair("EE",&eeloops_));
	loopsptr_.insert(std::make_pair("NH",&nhloops_));
	loopsptr_.insert(std::make_pair("NE",&neloops_));
	loopsptr_.insert(std::make_pair("HC",&hcloops_));
	loopsptr_.insert(std::make_pair("EC",&ecloops_));
}

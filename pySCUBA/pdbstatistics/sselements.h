/*
 * sselements.h
 *
 *  Created on: 2016年5月17日
 *      Author: hyliu
 */

#ifndef SSELEMENTS_H_
#define SSELEMENTS_H_
#include <backbone/backbonesite.h>
#include <map>
namespace NSPpdbstatistics {

class SSElements {
public:
//	enum LOOPS{HH=0,HE=1,EH=2,EE=3,NH=4,NE=5,HC=6,EC=7};
	SSElements(std::vector<NSPproteinrep::BackBoneSite> *sites);
	const std::vector<long> & helices(int l) const {
		if(helices_.find(l) != helices_.end())	return helices_.at(l);
		return empty;
	}
	const std::vector<long> & strands(int l) const {
		if(strands_.find(l) != strands_.end()) return strands_.at(l);
		return empty;
	}
	std::vector<std::vector<long>> & hhloops() {return hhloops_;}
	std::vector<std::vector<long>> & eeloops() {return eeloops_;}
	std::vector<std::vector<long>> & heloops() {return heloops_;}
	std::vector<std::vector<long>> & ehloops() {return ehloops_;}
	std::vector<std::vector<long>> & nhloops() {return nhloops_;}
	std::vector<std::vector<long>> & neloops() {return neloops_;}
	std::vector<std::vector<long>> & hcloops() {return hcloops_;}
	std::vector<std::vector<long>> & ecloops() {return ecloops_;}
	std::map<std::string,std::vector<std::vector<long>>*> & loopsptr() {return loopsptr_;}
	const std::vector<std::vector<long>> & hhloops() const {return hhloops_;}
	const std::vector<std::vector<long>> & eeloops() const {return eeloops_;}
	const std::vector<std::vector<long>> & heloops() const {return heloops_;}
	const std::vector<std::vector<long>> & ehloops() const {return ehloops_;}
	const std::vector<std::vector<long>> & nhloops() const {return nhloops_;}
	const std::vector<std::vector<long>> & neloops() const {return neloops_;}
	const std::vector<std::vector<long>> & hcloops() const {return hcloops_;}
	const std::vector<std::vector<long>> & ecloops() const {return ecloops_;}
	const std::map<std::string,std::vector<std::vector<long>>*> & loopsptr() const {return loopsptr_;}
	std::vector<NSPproteinrep::BackBoneSite> & sites(){return *sites_;}
private:
	std::vector<long> empty;
	std::vector<NSPproteinrep::BackBoneSite> *sites_;
	std::map<int,std::vector<long>>helices_;
	std::map<int,std::vector<long>> strands_;
	std::vector<std::vector<long>> hhloops_;
	std::vector<std::vector<long>> heloops_;
	std::vector<std::vector<long>> eeloops_;
	std::vector<std::vector<long>> ehloops_;
	std::vector<std::vector<long>> nhloops_; //nterminal pre helix loop
	std::vector<std::vector<long>> neloops_;
	std::vector<std::vector<long>> hcloops_; //cterminal after helix
	std::vector<std::vector<long>> ecloops_;
	std::map<std::string,std::vector<std::vector<long>>*> loopsptr_;
};
}



#endif /* SSELEMENTS_H_ */

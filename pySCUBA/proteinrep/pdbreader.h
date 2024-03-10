/*
 * readpdb.h
 *
 *  Created on: 2016年11月16日
 *      Author: hyliu
 */

#ifndef PROTEINREP_PDBREADER_H_
#define PROTEINREP_PDBREADER_H_
#include "dstl/mapkeyint.h"
#include "proteinrep/pdbrecord.h"
#include <cassert>
#include <map>
#include <memory>
namespace NSPproteinrep {
class MapPdbKeyInt;
class PdbReader {
public:
	typedef std::pair<int,char> ResKeyType; //residue number + insertion code
	typedef std::map<ResKeyType,std::vector<PdbRecord>> ResMapType;
	typedef std::map<char,ResMapType> RecordsType;
		// atom records organized by chain and residue
	static ResKeyType reskey(const std::string &residstr){
		int rid;
		char insertionid{' '};
		if(std::isdigit(residstr.back())){
			rid=std::stoi(residstr);
		} else{
			insertionid=residstr.back();
			rid=std::stoi(residstr.substr(0,residstr.size()-1));
		}
		return std::make_pair(rid,insertionid);
	}
	void readpdb(const std::string & filename);
	void readpdb(std::vector<std::string> & lines);
    RecordsType & records() {return records_;}
    const std::vector<std::string> &getaminoacidsequence(char chainid) const{
    	return aminoacidsequences_.at(chainid);
    }
    std::string chainids() const {
    	std::string res;
    	for(auto & c:records_) res.push_back(c.first);
    	return res;
    }
    const RecordsType & records() const {return records_;}
    void addRecord(const PdbRecord & record);
    std::shared_ptr<MapPdbKeyInt> mappdbkeyint() const {return mappdbkeyint_;}
private:
	RecordsType records_;
	std::map<char,std::vector<std::string>> aminoacidsequences_;
	std::shared_ptr<MapPdbKeyInt> mappdbkeyint_;
	 std::map<char,std::map<int, std::vector<char>>> multiinsertionids_; // occurred  insertionids at chain and positions
	 std::vector<char> remappedinsertionids{'a','b','c','d','e','f','g'};
};

class MapPdbKeyInt {
public:
	MapPdbKeyInt(const typename PdbReader::RecordsType & records,
			const std::map<char,std::map<int, std::vector<char>>> & multiinsertionids,
			const std::vector<char> & remappedids);
	char pdbChainID(unsigned int chainnumber=0)
		{return mapchainidint_.key(chainnumber);}
	typename PdbReader::ResKeyType pdbResKey(unsigned int posinumber,unsigned int chainnumber=0) {
		return mapreskeyint_.at(chainnumber).key(posinumber);  //remapped insertionid, not original insertion id
	}
	int pdbResID(unsigned int posinumber, unsigned int chainnumber=0){
		return pdbResKey(posinumber, chainnumber).first;
	}
	unsigned int chainNumber(char pdbchainid){
		return mapchainidint_.keynumber(pdbchainid);
	}
	unsigned int posiNumber(const PdbReader::ResKeyType &reskey_or, char pdbchainid=' ') {
		char mappedid=reskey_or.second;
			if(multiinsertionids_.find(pdbchainid) != multiinsertionids_.end()){
				if(multiinsertionids_[pdbchainid].find(reskey_or.first) != multiinsertionids_[pdbchainid].end()){
					auto iter =remappedinsertionids.begin();
					for(char t:multiinsertionids_[pdbchainid][reskey_or.first]){
	                     if(reskey_or.second==t) {
	                    	 mappedid=*iter;
	                    	 break;
	                     }
	                     ++iter;
					}
				}
			}
		return mapreskeyint_.at(chainNumber(pdbchainid)).keynumber(std::make_pair(reskey_or.first,mappedid));
	}
   const NSPdstl::MapKeyInt<char> & mapchainidint() const {return mapchainidint_;}
/*	const NSPdstl::MapKeyInt<typename PdbReader::ResKeyType> mapreskeyint(unsigned int chainNumber) const {
		return mapreskeyint_.at(chainNumber);
	}*/

private:
	NSPdstl::MapKeyInt<char> mapchainidint_;
	std::vector<NSPdstl::MapKeyInt<typename PdbReader::ResKeyType>> mapreskeyint_;
	 std::map<char,std::map<int, std::vector<char>>> multiinsertionids_; // occurred  insertionids at chain and positions
	 std::vector<char> remappedinsertionids;
};

struct PDBModelID {
	std::string pdbid{""};
	int biounit{1};
	int model{0};
};


std::string extractpdbmodel(const PDBModelID & mid);
std::string extractpdbmodel(const std::string & filename,const PDBModelID &mid);
std::shared_ptr<PdbReader> readpdbmodel(const std::string & pdbid,int biounit=1,
		int modelno=0);
}

#endif /* PROTEINREP_PDBREADER_H_ */

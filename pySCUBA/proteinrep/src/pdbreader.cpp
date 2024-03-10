/*
 * pdbreader.cpp
 *
 *  Created on: 2016年11月16日
 *      Author: hyliu
 */

#include "proteinrep/pdbreader.h"
#include "dataio/inputlines.h"
#include "dataio/datapaths.h"
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace NSPproteinrep;
void PdbReader::readpdb(const std::string & filename) {
	NSPdataio::TextLines lines;
	lines.init(filename);
	readpdb(lines.lines());
}

void PdbReader::readpdb(std::vector<std::string> & lines) {
	for (auto & line : lines) {
		if (line.substr(0, 6) != "ATOM  " && line.substr(0, 6) != "HETATM")
			continue;
		addRecord(PdbRecord(line));
	}
	mappdbkeyint_=std::shared_ptr<MapPdbKeyInt>(
			new MapPdbKeyInt(records_,multiinsertionids_,remappedinsertionids));
	for( auto &c:records_) {
		aminoacidsequences_.insert(std::make_pair(c.first,std::vector<std::string>()));
		std::vector<std::string>&seq=aminoacidsequences_.at(c.first);
		int posi=0;
		for(auto &r:c.second){
//			if(r.second[0].label != "ATOM") continue;
			assert(posi++ == mappdbkeyint_->posiNumber({r.second[0].residueid,r.second[0].insertionid},c.first));
			seq.push_back(r.second[0].residuename);
		}
	}
}

void PdbReader::addRecord(const PdbRecord & record) {
	char chainid = record.chainid;
	if (records_.find(chainid) == records_.end())
		records_.insert(std::make_pair(chainid, ResMapType()));
	ResMapType & resmap = records_.at(chainid);
	char mappedid=record.insertionid;
	bool noinsertion=record.insertionid==' ';
	if(noinsertion) {
		if(this->multiinsertionids_.find(record.chainid)!= this->multiinsertionids_.end()){
			noinsertion=multiinsertionids_[chainid].find(record.residueid)==multiinsertionids_[chainid].end();
		}
	}
	if(!noinsertion){
		if(this->multiinsertionids_.find(record.chainid)== this->multiinsertionids_.end()){
			multiinsertionids_[chainid]=std::map<int,std::vector<char>>();
		}
		if(multiinsertionids_[chainid].find(record.residueid)==multiinsertionids_[chainid].end()){
			multiinsertionids_[chainid][record.residueid]=std::vector<char>();
		}
	    bool notadded=true;
		for(auto iid:multiinsertionids_[chainid][record.residueid]) {
			if (iid== record.insertionid) notadded=false;
		}
		if(notadded) multiinsertionids_[chainid][record.residueid].push_back(record.insertionid);
		int iid;
		for(int i=0;i<multiinsertionids_[chainid][record.residueid].size();++i){
			if(record.insertionid==multiinsertionids_[chainid][record.residueid][i]) iid=i;
		}
		mappedid=this->remappedinsertionids[iid];
//		std::cout <<"record after insertion " <<record.toString()<< "mapped code : " << mappedid<<std::endl;
	}
	ResKeyType reskey = std::make_pair(record.residueid, mappedid);
	if (resmap.find(reskey) == resmap.end())
		resmap.insert(std::make_pair(reskey, std::vector<PdbRecord>()));
	resmap.at(reskey).push_back(record);
}

MapPdbKeyInt::MapPdbKeyInt(const typename PdbReader::RecordsType & records,
		const std::map<char,std::map<int, std::vector<char>>> & multiinsertionids,
		const std::vector<char> & remappedids) {
	multiinsertionids_=multiinsertionids;
	remappedinsertionids=remappedids;
	mapchainidint_.init(records);
	mapreskeyint_.resize(records.size(),
			NSPdstl::MapKeyInt<typename PdbReader::ResKeyType>());
	for (auto & c : records)
		mapreskeyint_.at(mapchainidint_.keynumber(c.first)).init(c.second);
}
std::string pdbfilename(const PDBModelID &mid) {
	std::string res;
	res="pdb"+mid.pdbid+".ent";
/*	res.resize(mid.pdbid.size());
	std::transform(mid.pdbid.begin(), mid.pdbid.end(), res.begin(), ::tolower);
	res += ".pdb" + std::to_string(mid.biounit);*/
	return res;
}
std::string modelfilename(const PDBModelID &mid) {
	std::string res;
	res = mid.pdbid + "_model" + std::to_string(mid.model) + ".pdb"+std::to_string(mid.biounit);
	return res;
}
static std::vector<std::string>  extractmodellines(const std::string &filename,int modelno){
	std::vector<std::string> res;
	NSPdataio::TextLines lines;
	lines.init(filename);
	bool inmodel = false;
	if (modelno == 0 ||modelno ==1)
		inmodel = true;
	int natoms = 0;
	for (auto &line : lines) {
		if (line.substr(0, 6) != "ATOM  " && line.substr(0, 6) != "HETATM"
				&& line.substr(0, 6) != "ENDMDL"
				&& line.substr(0, 5) != "MODEL") {
			res.push_back(line);
			continue;
		}
		if (line.substr(0, 5) == "MODEL") {
			inmodel = false;
			int num = std::stoi(line.substr(6));
			if (num == modelno)
				inmodel = true;
			continue;
		}
		if (line.substr(0, 6) == "ENDMDL") {
			inmodel = false;
			continue;
		}
		if (inmodel) {
			res.push_back( line);
			++natoms;
		}
	}
	return res;
}
static std::string fndordldpdbfile(const PDBModelID &mid){
	std::string pdbpath = NSPdataio::downloadedpdbpath();
	std::string filename = pdbfilename(mid);
	std::string pfname=pdbpath+filename;
	return pfname;
}
std::string NSPproteinrep::extractpdbmodel(const PDBModelID & mid) {
	std::string filename=fndordldpdbfile(mid);
	return extractpdbmodel(filename, mid);
}

std::string NSPproteinrep::extractpdbmodel(const std::string & filename,
		const PDBModelID &mid) {
	std::vector<std::string> lines=extractmodellines(filename,mid.model);
	if(lines.empty()) return std::string();
	std::string mfile = modelfilename(mid);
	std::ofstream ofs;
	ofs.open(mfile.c_str());
	for(auto &l:lines) ofs <<l<<std::endl;
	return mfile;
	ofs.close();
}
std::shared_ptr<PdbReader> NSPproteinrep::readpdbmodel(const std::string & pdbid,int biounit,
		int modelno){
	PDBModelID mid;
	mid.pdbid=pdbid;
	mid.biounit=biounit;
	mid.model=modelno;
	std::string pdbfile=fndordldpdbfile(mid);
	if(pdbfile.empty()) return nullptr;
	std::vector<std::string> lines=extractmodellines(pdbfile,modelno);
	std::string cmd = "gzip " + pdbfile;
	FILE *pp=popen(cmd.c_str(), "w");
	std::fclose(pp);
	std::shared_ptr<PdbReader> reader=std::shared_ptr<PdbReader> (new PdbReader());
	reader->readpdb(lines);
	return reader;
}
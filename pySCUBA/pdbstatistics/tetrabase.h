/*
 * rminsef.h
 *
 *  Created on: 2016年6月8日
 *      Author: hyliu
 */

#ifndef RMINSEF_H_
#define RMINSEF_H_
#include <backbone/backbonesite.h>
#include <map>

namespace NSPpdbstatistics {
class BBD2Matrix {
public:
	enum {N=0,C=1,CB=2,CX=3};
	BBD2Matrix(NSPproteinrep::BackBoneSite &s1, NSPproteinrep::BackBoneSite &s2);
	BBD2Matrix(NSPproteinrep::BackBoneSite &s1, NSPproteinrep::BackBoneSite &s2, NSPproteinrep::BBatomdis2 &atomdis2);
	int min_index() const{ return min_index_;}
	int max_index() const{ return max_index_;}
	std::pair<int,int> minmaxindices() const {return std::make_pair(min_index_,max_index_);}
	std::vector<double> &matrix(){return dist2_;}
	double phi() {return phi_;}
private:
	std::vector<double> dist2_;
	int min_index_;
	int max_index_;
	double phi_;
	std::pair<int,int>indextopair(int idx);
};
class EneTables {
public:
	struct Table {
		std::vector<double> r;
		std::vector<double> ene;
	};
	typedef std::pair<int,int> OrienType;
	typedef std::string SSType;
	typedef std::pair<SSType,OrienType> EneType;

	static double energy(const Table &tab, double r);
	static double subenergy(const Table &rtab,const std::vector<Table> &ptables,double r, double p);
/*	double energy(EneType & t, double r) const {
		if(tables_.find(t) == tables_.end()) return 10.0;
		return energy(tables_.at(t),r);
	}
	*/
	double energy(EneType & t, double r,double p) const {
		if(tables_.find(t) == tables_.end()) return 10.0;
		double etot=energy(tables_.at(t),r);
		if(phiene_) etot +=subenergy(tables_.at(t),subtables_.at(t),r,p);
		return etot;
	}
	std::map<EneType,Table> & tables() {return tables_;}
	const std::map<EneType,Table> & tables() const {return tables_;}
//	void read(const std::string & filename);
	void readwithsubtables(const std::string &filename);
private:
	static int getindex(const std::vector<double> &crds,double c) {
		return intervale(crds,c,0,crds.size()-1);
	}
	static int intervale(const std::vector<double> &crds,double c, int lower, int upper);
	std::map<EneType,Table> tables_;
	std::map<EneType,std::vector<Table>> subtables_;
	bool phiene_{false};
};

class TetraBASE {
public:
	void init(const std::string & filename) {
			enetables_.readwithsubtables(filename);
	}
	static typename EneTables::EneType enetype(char s1, char s2, const BBD2Matrix & d2m);
	double onebody(NSPproteinrep::BackBoneSite &s) {return 0.0;}
	double twobody(NSPproteinrep::BackBoneSite &s1, NSPproteinrep::BackBoneSite &s2);
	double twobody(NSPproteinrep::BackBoneSite &s1, NSPproteinrep::BackBoneSite &s2,NSPproteinrep::BBatomdis2 &atomdis2);
	EneTables & enetables() {return enetables_;}
	const EneTables & enetables() const {return enetables_;}
private:
	EneTables enetables_;
};
}

#endif /* RMINSEF_H_ */

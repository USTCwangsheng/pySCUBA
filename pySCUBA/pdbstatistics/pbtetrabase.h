/*
 * pbtetrabase.h
 *
 *  Created on: 2017年8月4日
 *      Author: hyliu
 */

#ifndef PDBSTATISTICS_PBTETRABASE_H_
#define PDBSTATISTICS_PBTETRABASE_H_
#include "backbone/backbonesite.h"
#include "backbone/hbonded.h"
#include <map>

namespace NSPpdbstatistics {


class TetraGeom {
public:
	typedef std::pair<int, int> OrientType;
	typedef std::pair<char, char> PbTypes;
	typedef std::pair<PbTypes, OrientType> KeyType;
	void init(const std::vector<double> &dist2);
	TetraGeom(const NSPproteinrep::BackBoneSite &s1,
			const NSPproteinrep::BackBoneSite &s2,
			const std::vector<double> & dist2) {
		char sscode1=s1.sscode;
		if(sscode1=='H') sscode1='m';
		if(sscode1=='E') sscode1='d';
		char sscode2=s2.sscode;
		if(sscode2=='H') sscode2='m';
		if(sscode2=='E') sscode2='d';
		pbtypes_ = std::make_pair(sscode1, sscode2);
		init(dist2);
		if(dist2[1]<22.4 || dist2[4]<22.4){
			hbond_=(NSPproteinrep::hbonded(s1,s2,&hbg_)!=0);
		}
	}
	TetraGeom(const NSPproteinrep::BackBoneSite &s1,
			const NSPproteinrep::BackBoneSite &s2) {
		std::vector<double> dist2;
		char sscode1=s1.sscode;
		if(sscode1=='H') sscode1='m';
		if(sscode1=='E') sscode1='d';
		char sscode2=s2.sscode;
		if(sscode2=='H') sscode2='m';
		if(sscode2=='E') sscode2='d';
		pbtypes_ = std::make_pair(sscode1, sscode2);;
		dist2matrix(s1, s2, &dist2);
		init(dist2);
		if(dist2[1]<22.4 ||dist2[4]<22.4){
			hbond_=(NSPproteinrep::hbonded(s1,s2,&hbg_)!=0);
		}
	}
	KeyType getkey() const {
		return std::make_pair(pbtypes_, orient_);
	}
	OrientType orient() const {
		return orient_;
	}
	PbTypes pbtypes() const {
		return pbtypes_;
	}
	double rmin() const {
		return rmin_;
	}
	bool hbond() const { return hbond_;}
	const NSPproteinrep::HBondGeometry &hbbondgeometry() const {return hbg_;}
private:
	OrientType orient_;
	PbTypes pbtypes_;
	double rmin_;
	bool hbond_{false};
	NSPproteinrep::HBondGeometry hbg_;
	void dist2matrix(const NSPproteinrep::BackBoneSite &s1,
			const NSPproteinrep::BackBoneSite &s2, std::vector<double> *dist2);
};
template<typename ENETYPE>
class ETables {
public:
	typedef ENETYPE EneType;
	struct Table {
		std::vector<double> r;
		std::vector<double> ene;
	};
	static double energy(const Table &tab, double r){
		const std::vector<double> & crds=tab.r;
		const std::vector<double> & values=tab.ene;
//		if(r>crds.back()) r=crds.back();
		if(r>crds.back()) return 0.0;
		int idx=getindex(crds,r);
		double alpha=(r-crds[idx])/(crds[idx+1]-crds[idx]);
		return values[idx]*alpha +values[idx+1]*(1-alpha);
	}
	double energy(const EneType & etype, double r) const {
		if (tables_.find(etype) == tables_.end())
			return 100.0;
		double etot = energy(tables_.at(etype), r);
		return etot;
	}
	void inserttable(EneType &etype, const std::vector<double> & r,
			const std::vector<double> &ene){
		tables_.insert(std::make_pair(etype,Table()));
		tables_[etype].r=r;
		tables_[etype].ene=ene;
	}
	bool empty() const {return tables_.empty();}
private:
	static int getindex(const std::vector<double> &crds, double c) {
		return intervale(crds, c, 0, crds.size() - 1);
	}
	static int intervale(const std::vector<double> &crds, double c, int lower,
			int upper) {
		if (upper <= lower + 1)
			return lower;
		int middle = lower + (upper - lower) / 2;
		if (c < crds[middle])
			upper = middle;
		else
			lower = middle;
		return intervale(crds, c, lower, upper);
	}
	std::map<EneType, Table> tables_;
};
typedef std::pair<TetraGeom::PbTypes, TetraGeom::OrientType> TetraPairType;
typedef ETables<TetraPairType> TetraPairETables;
bool readetables(const std::string & filename, TetraPairETables *);
bool readetables(const std::string & filename, TetraPairETables * Ermin,
		TetraPairETables *Ehb_cor, TetraPairETables *eno_hb_cor);
TetraPairETables & gettetrapairetables(const std::string &filename="");
TetraPairETables & gettetrapairetables(bool hbond, const std::string &filename="");
double tetrapairenergy(const NSPproteinrep::BackBoneSite &s1, const NSPproteinrep::BackBoneSite &s2);
}

#endif /* PDBSTATISTICS_PBTETRABASE_H_ */

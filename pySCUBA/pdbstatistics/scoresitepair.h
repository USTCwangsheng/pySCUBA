/*
 * scoresitepair.h
 *
 *  Created on: 2017年12月5日
 *      Author: hyliu
 */

#ifndef PDBSTATISTICS_SCORESITEPAIR_H_
#define PDBSTATISTICS_SCORESITEPAIR_H_
#define USE_CB

#include "backbone/backbonesite.h"
#include "dstl/vectortree.h"
#include "dstl/randomengine.h"

namespace NSPpdbstatistics{
typedef  std::pair<std::vector<NSPproteinrep::BackBoneSite>::const_iterator,
				   std::vector<NSPproteinrep::BackBoneSite>::const_iterator> SiteItPair;
typedef std::pair<NSPproteinrep::BackBoneSite,NSPproteinrep::BackBoneSite> SitePair;
class QueryPair{
public:
	QueryPair(){;}
	QueryPair(const SiteItPair &conf);
	QueryPair swapped() const{
		QueryPair result;
		result.ls1()=ls2_;
		result.ls2()=ls1_;
		result.crd2()=crd1_;
		result.crd1()=crd2_;
		result.alpha()=alpha_;
		result.ncut1()=ncut2_;
		result.ncut2()=ncut1_;
		return result;
	}
	std::vector<double> dmvector() const;
	double dca2() const{return (crd1_[1]-crd2_[1]).squarednorm();}
	std::vector<double> &ls1(){return ls1_;}
	std::vector<double> &ls2(){return ls2_;}
	std::vector<NSPgeometry::XYZ> &crd1(){return crd1_;}
	std::vector<NSPgeometry::XYZ> &crd2(){return crd2_;}
	double & alpha() {return alpha_;}
	int &ncut1(){return ncut1_;}
	int &ncut2(){return ncut2_;}
	std::vector<std::vector<double>> & scores() {return scores_;}
	std::vector<double> summscore(const std::vector<double> &tmpltsizes,
			const std::vector<double> &refsizes,double refblank=0.0001);
	void saveconf(std::ofstream &ofs);
private:
	std::vector<double> ls1_; //local torsions of site1
	std::vector<double> ls2_; //local torsions of site2
	std::vector<NSPgeometry::XYZ> crd1_; //atomic coordinates of site1
	std::vector<NSPgeometry::XYZ> crd2_; //atomic coordinates of site2
	double alpha_{-1.0}; //parameter for searching distance matrix neighbors;
	int ncut1_{-1}; //parameter for searching local structure neighbors
	int ncut2_{-1}; //parameter for searching local structure neighbors;
	std::vector<std::vector<double>> scores_; // scores computed using different template groups
};
class ScoreSitePair{
public:
	ScoreSitePair(double dmin,double dmax):dcamin_(dmin),dcamax_(dmax){;}
	void builddmtree(const std::vector<NSPproteinrep::BackBoneSite> & sites,bool countpaironly=false);
	void buildlstree(const std::vector<NSPproteinrep::BackBoneSite> & sites,double pkeep);
	void buildrefdmtree(int ntimes);
	void score(QueryPair &qp,bool silent=false);
	const std::vector<SiteItPair> & tmpltsitepairs() const {return sips_;}
	double neighborsum(QueryPair &qp);
//	double neighborsum(const SiteItPair &conf, double alpha,int nlscut1,int nlscut2);
	double refdmneighborsum(QueryPair &qp);
	double reflsneighborsum(QueryPair &qp);

	static std::pair<double,double> coverrange(double dcamin,double dcamax){
		double dmin=dcamin+0.6;
		if(dcamin<0.01) dmin=0.0;
		return std::make_pair(dmin,dcamax-0.6);
	}

	bool covered(double dca) const {
		auto range=coverrange(dcamin_,dcamax_);
		return (dca>=range.first &&dca<range.second);}

	bool covered(const SiteItPair &conf) const {
		double rca=sqrt((conf.first->cacrd()-conf.second->cacrd()).squarednorm());
		return covered(rca);
	}
	long tmpltsize() const {return dms_->size();}
	long refsize() const {return refdms_->size();}
private:
	double dcamin_;
	double dcamax_;
	domaintree::VectorTree dmtree_;
	domaintree::VectorTree refdmtree_;
	domaintree::VectorTree lstree_;
	std::vector<SiteItPair> sips_;
	std::shared_ptr<std::vector<std::vector<double>>> dms_;
	std::shared_ptr<std::vector<std::vector<double>>> refdms_;
	std::shared_ptr<std::vector<std::vector<double>>> refls_;
	static double dmmatchscore(const std::vector<double> & dm, const std::vector<double> &dmtmplt,double alpha);
	static double lsmatchscore(const std::vector<double> &ls, const std::vector<double> &lstmplt);

	static std::vector<double> dmvector(const NSPproteinrep::BackBoneSite &s1,
			const NSPproteinrep::BackBoneSite &s2);
};

std::vector<NSPproteinrep::BackBoneSite> makerandompairconf(const SiteItPair &conf,
		int ntime,double dcamin,double dcamax,NSPdstl::RandomEngine<> & reng);

inline std::vector<NSPproteinrep::BackBoneSite> makerandompairconf(const SiteItPair &conf,
		int ntime,double dcamin,double dcamax){
	auto & reng=NSPdstl::RandomEngine<>::getinstance();
	return makerandompairconf(conf,ntime,dcamin,dcamax,reng);
}

std::vector<long> splitsites(const std::vector<NSPproteinrep::BackBoneSite> &sites,
		int ngroups);

void drawquerypairs(const std::vector<NSPproteinrep::BackBoneSite> &sites,
		double dcamin,
		double dcamax,
		int nskip, int nsample,int nrandomtime,
		std::vector<QueryPair> *qp_native, std::vector<QueryPair> *qp_random);
}


#endif /* PDBSTATISTICS_SCORESITEPAIR_H_ */

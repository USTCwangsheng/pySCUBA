/*
 * scoresitepair.cpp
 *
 *  Created on: 2017年12月5日
 *      Author: hyliu
 */
#include "pdbstatistics/scoresitepair.h"
#include "dstl/randomengine.h"

using namespace NSPproteinrep;
using namespace NSPpdbstatistics;
using namespace domaintree;
using namespace NSPgeometry;
#include <cmath>
const double ALPHA0 = 0.01;
//double alpha;
std::vector<double> lscut { 40, 20, 20, 40 };
static std::vector<std::pair<double, double>> dcuts(
		const std::vector<double> & dm, double alpha) {
	std::vector<std::pair<double, double>> res;
	for (auto d : dm) {
		double dcut1 = alpha * (d * d);
		double dcut2 = 2.0 * dcut1;
		res.push_back(std::make_pair(dcut1, dcut2));
	}
	return res;
}
static double bound2(const std::vector<double> &dm, double alpha) {
	double res = 0.0;
	for (auto d : dm) {
		double dcut2 = 2.0 * alpha * (d * d);
		res += dcut2 * dcut2;
	}
	return res;
}
/*
static bool steric_clash(const std::vector<double> dm) {
	static const std::vector<double> clashdistance { 3.0, 3.0, 3.0, 2.5, 3.0,
			3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 2.5, 3.0, 3.0, 3.0 };
	for (int i = 0; i < dm.size(); ++i) {
		if (dm[i] < clashdistance[i])
			return true;
	}
	return false;
}*/

static std::vector<double> getls(
		const std::vector<NSPproteinrep::BackBoneSite>::const_iterator & it) {
	std::vector<double> ls;
	ls.push_back((it - 1)->psi());
	ls.push_back(it->phi());
	ls.push_back(it->psi());
	ls.push_back((it+1)->phi());
	for (auto &a : ls) {
		if (a >= 360.0)
			a -= 360.0;
		if (a <= 0)
			a += 360.0;
	}
	return ls;
}
QueryPair::QueryPair(const SiteItPair &conf) {
	static const std::vector<int> crdidx { BackBoneSite::NCRD,
			BackBoneSite::CACRD, BackBoneSite::CCRD, BackBoneSite::OCRD };
	ls1_ = getls(conf.first);
	ls2_ = getls(conf.second);
	for (int i = 0; i < 4; ++i) {
		crd1_.push_back(conf.first->getcrd(crdidx[i]));
		crd2_.push_back(conf.second->getcrd(crdidx[i]));
	}
#ifdef USE_CB
	crd1_.push_back(conf.first->cbcrd());
	crd2_.push_back(conf.second->cbcrd());
#endif
}

std::vector<double> QueryPair::dmvector() const {
	std::vector<double> dm;
	for (auto &p1 : crd1_) {
		for (auto &p2 : crd2_)
			dm.push_back(sqrt((p1 - p2).squarednorm()));
	}
	return dm;
}
std::vector<double> QueryPair::summscore(const std::vector<double> &tmpltsizes,
		const std::vector<double> &refsizes, double refblank) {
	double stot = 0.0;
	double srtot = 0.0;
	double tstot = 0.0;
	double rstot = 0.0;
	std::vector<double> res(1 + 2 * scores_.size());
	for (int i = 0; i < scores_.size(); ++i) {
		res[2 * i + 1] = scores_[i][0];
		res[2 * i + 2] = scores_[i][1];
		stot += scores_[i][0];
		srtot += scores_[i][1];
		tstot += tmpltsizes[i];
		rstot += refsizes[i];
	}
	double mins = refblank / rstot;
	res[0] = (stot / tstot + mins) / (srtot / rstot + mins);
	return res;
}
void QueryPair::saveconf(std::ofstream &ofs){
	for(auto t:ls1_) ofs <<" "<<t;
	ofs<<std::endl;
	for(auto x:crd1_) ofs<<x.toString()<<std::endl;
	for(auto t:ls2_) ofs <<" "<<t;
	ofs<<std::endl;
	for(auto x:crd2_) ofs<<x.toString()<<std::endl;
}
static double matchls(const std::vector<double> &ls1,
		const std::vector<double> &ls2, int ncut) {
	double res = 1.0;
	std::vector<double> lscut1(4), lscut2(4);
	for (int i = 0; i < 4; ++i) {
		lscut1[i] = pow(1.5, ncut) * lscut[i];
		lscut2[i] = 1.5 * lscut1[i];
	}
	for (int i = 0; i < 4; ++i) {
		double diff = ls1[i] - ls2[i];
		if (diff > 180.0)
			diff -= 360.0;
		if (diff < -180.0)
			diff += 360.0;
		diff = fabs(diff);
		if (diff <= lscut1[i])
			continue;
		if (diff > lscut2[i]) {
			res = 0.0;
			break;
		}
		double x = (diff - lscut2[i]) / (lscut1[i] - lscut2[i]);
		res *= (x * x);
	}
	return res;
}
static double lsbound2(int ncut) {
	std::vector<double> lscut2(4);
	for (int i = 0; i < 4; ++i) {
		lscut2[i] = 1.5 * pow(1.5, ncut) * lscut[i];
	}
	double res = 0.0;
	for (auto t : lscut2)
		res += t * t;
	return res;
}
std::vector<double> ScoreSitePair::dmvector(const BackBoneSite &s1,
		const BackBoneSite &s2) {
	std::vector<double> dm;
	std::vector<XYZ> crd1;
	s1.getcrd(crd1);
	std::vector<XYZ> crd2;
	s2.getcrd(crd2);
#ifdef USE_CB
	crd1.push_back(s1.cbcrd());
	crd2.push_back(s2.cbcrd());
#endif
//	int idx = 0;
	for (auto &p1 : crd1) {
		for (auto &p2 : crd2)
			dm.push_back(sqrt((p1 - p2).squarednorm()));
	}
	return dm;
}
void ScoreSitePair::builddmtree(
		const std::vector<NSPproteinrep::BackBoneSite> & sites,
		bool countpaironly) {
	int minsep = 5;
	auto start2 = sites.begin() + 1;
	long npairs = 0;
	if (!countpaironly)
		dms_ =
				std::shared_ptr < std::vector<std::vector<double>>>(new std::vector<
						std::vector<double>>());
	for (auto iter1 = start2 + minsep; iter1 != sites.end() - 1; ++iter1) {
		if (chainstartsite(iter1)) {
//			std::cout<<"Number of chains: " <<++nchain <<" Number of Pairs: " <<paircount <<std::endl;
			start2 = iter1 + 1;
		}

		if (iter1 - start2 < minsep)
			continue;
		if (!fragstartsite(iter1 - 1, sites.end(), 3))
			continue;
		if (iter1->resname=="GLY") continue;
		int sslength = 0;
		if (iter1->sscodechar() != 'C') {
			while (true) {
				if (iter1 - sslength == start2)
					break;
				if ((iter1 - (sslength + 1))->sscodechar()
						== iter1->sscodechar()) {
					++sslength;
					continue;
				} else {
					break;
				}
			}
		}
		BackBoneSite s1 = *iter1;
		int iter2end = sslength > minsep ? sslength : minsep;
		for (auto iter2 = start2; iter2 != iter1 - iter2end; ++iter2) {
			if (!fragstartsite(iter2 - 1, sites.end(), 3))
				continue;
			if(iter2->resname=="GLY") continue;
			if(s1.resname=="ALA" && iter2->resname=="ALA") continue;
			BackBoneSite s2 = *iter2;
			double dca2 = (s2.cacrd() - s1.cacrd()).squarednorm();
			if (dca2 < dcamin_ * dcamin_ || dca2 >= dcamax_ * dcamax_)
				continue;
			std::vector<double> dm = dmvector(s2, s1);
			if (dm.empty())
				continue;
			if (!countpaironly) {
				sips_.push_back(std::make_pair(iter2, iter1));
				dms_->push_back(dm);
				if (npairs == 0) {
					dmtree_.init(dms_, 0.1, 0.0, 20.0);
				}
				dmtree_.insertpoint(npairs);
			}
			++npairs;
		}
	}
	std::cout << "Number of site pairs: " << npairs << std::endl;
}
void ScoreSitePair::buildlstree(
		const std::vector<NSPproteinrep::BackBoneSite> & sites, double pkeep) {
	int minsep = 5;
	auto start2 = sites.begin() + 1;
	long nls = 0;
	refls_ =
			std::shared_ptr < std::vector<std::vector<double>>>(new std::vector<
					std::vector<double>>());
	for (auto iter1 = start2 + minsep; iter1 != sites.end() - 1; ++iter1) {
		if (chainstartsite(iter1)) {
//			std::cout<<"Number of chains: " <<++nchain <<" Number of Pairs: " <<paircount <<std::endl;
			start2 = iter1 + 1;
		}

		if (iter1 - start2 < minsep)
			continue;
		if (!fragstartsite(iter1 - 1, sites.end(), 3))
			continue;
		double rn = NSPdstl::RandomEngine<>::getinstance().realrng(0.0, 1.0)();
		if (rn > pkeep)
			continue;
		std::vector<double> ls = getls(iter1);
		refls_->push_back(ls);
		if (nls == 0) {
			lstree_.init(refls_, 1.0, 0.0, 360.0);
		}
		lstree_.insertpoint(nls);
		++nls;
	}
}
double ScoreSitePair::dmmatchscore(const std::vector<double> &dm,
		const std::vector<double> &tmpltdm, double alpha) {
	std::vector<std::pair<double, double>> dcut = dcuts(dm, alpha);
	double s = 1.0;
	for (int i = 0; i < dm.size(); ++i) {
		double d = fabs(dm[i] - tmpltdm[i]);
		if (d < dcut[i].first)
			continue;
		if (d >= dcut[i].second) {
			s = 0.0;
			break;
		}
		double x = (d - dcut[i].second) / (dcut[i].first - dcut[i].second);
		s *= (x * x);
	}
	return s;
}

double ScoreSitePair::neighborsum(QueryPair &qp) {
	std::vector<double> dm = qp.dmvector();
	assert(qp.alpha() > 0.0);
	assert(qp.ncut1() >= 0);
	assert(qp.ncut2() >= 0);
	double bnd2 = bound2(dm, qp.alpha());
	D2Leaf<long, std::vector<std::vector<double>>, UsualCrd> d2leaf(dms_.get(),
			1000000, bnd2);
	dmtree_.gettree().findneighbor(dm, d2leaf, bnd2);
	std::vector<std::pair<long, double>> &neighbors =
			d2leaf.nnearest().neighbors();
	double s = 0.0;
//	double s=NMIN;
	for (auto &n : neighbors) {
		if (n.second < 0.000001)
			continue;  //ignore self
		const std::vector<double> & temp = dms_->at(n.first);
		double sdm = dmmatchscore(dm, temp, qp.alpha());
		std::vector<double> lstmp = getls(sips_.at(n.first).first);
		double sls1 = matchls(qp.ls1(), lstmp, qp.ncut1());
		lstmp = getls(sips_.at(n.first).second);
		double sls2 = matchls(qp.ls2(), lstmp, qp.ncut2());
		s += sdm * sls1 * sls2;
	}
	return s;
}
std::vector<BackBoneSite> NSPpdbstatistics::makerandompairconf(
		const SiteItPair &conf, int ntime, double dcamin, double dcamax,
		NSPdstl::RandomEngine<> & reng) {
	std::vector<BackBoneSite> res(3);
	std::copy(conf.first - 1, conf.first + 2, res.begin());
	std::vector<BackBoneSite> pep2(3);
	std::copy(conf.second - 1, conf.second + 2, pep2.begin());
	for (int i = 0; i < 3; ++i)
		res[i].translate(-1.0 * res[1].cacrd());
	NSPgeometry::XYZ trans;
	int nget = 0;
	while (nget < ntime) {
		double norm = -1.0;
		while (norm <= dcamin) {
			trans = NSPgeometry::XYZ(reng.realrng(0.0, 1.0), dcamax);
			norm = sqrt(trans.squarednorm());
		}
		trans = trans - pep2[1].cacrd();
		for (int i = 0; i < 3; ++i)
			pep2[i].translate(trans);
		NSPgeometry::XYZ axis(reng.realrng(), 1.0);
		double angle = reng.realrng(0.0, 1.0)() * 180.0;
		NSPgeometry::Rotation rot;
		rot.init(NSPgeometry::QuaternionCrd(axis, angle), pep2[1].cacrd());
		for (int i = 0; i < 3; ++i)
			pep2[i].rotate(rot);
		bool clashed = false;
		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j) {
				if (NSPproteinrep::atomsclashed(res[i], pep2[j])) {
					clashed = true;
					break;
				}
			}
			if (clashed)
				break;
		}
		if (!clashed) {
			for (int i = 0; i < 3; ++i) {
				res.push_back(pep2[i]);
			}
			++nget;
		}
	}
	return res;
}
void ScoreSitePair::buildrefdmtree(int ntimes) {
	refdms_ =
			std::shared_ptr < std::vector<std::vector<double>>>(new std::vector<
					std::vector<double>>());
	auto &reng = NSPdstl::RandomEngine<>::getinstance();
	int npairs = 0;
	for (auto & sip : sips_) {
		std::vector<BackBoneSite> randompairs = makerandompairconf(sip, ntimes,
				dcamin_, dcamax_);
		assert(randompairs.size() == 3 + 3 * ntimes);
		for (int j = 0; j < ntimes; ++j) {
			std::vector<double> dm = dmvector(randompairs[1],
					randompairs[3 * j + 4]);
			refdms_->push_back(dm);
			if (npairs == 0) {
				refdmtree_.init(refdms_, 0.1, 0.0, 20.0);
			}
			refdmtree_.insertpoint(npairs);
			++npairs;
		}
	}
}
double ScoreSitePair::refdmneighborsum(QueryPair &qp) {
	std::vector<double> dm = qp.dmvector();
	assert(qp.alpha() > 0.0);
	double bnd2 = bound2(dm, qp.alpha());
	D2Leaf<long, std::vector<std::vector<double>>, UsualCrd> d2leaf(
			refdms_.get(), 1000000, bnd2);
	refdmtree_.gettree().findneighbor(dm, d2leaf, bnd2);
	std::vector<std::pair<long, double>> &neighbors =
			d2leaf.nnearest().neighbors();
	double s = 0.0;
//	double s=NMIN;
	for (auto &n : neighbors) {
		if (n.second < 0.000001)
			continue;  //ignore self
		const std::vector<double> & temp = refdms_->at(n.first);
		s += dmmatchscore(dm, temp, qp.alpha());
	}
	return s;
}
double ScoreSitePair::reflsneighborsum(QueryPair &qp) {
	std::vector<double> &ls1 = qp.ls1();
	std::vector<double> &ls2 = qp.ls2();
	double s1;
	if (qp.ncut1() < 0) {
		qp.ncut1() = -1;
		do {
			++(qp.ncut1());
			double bnd = lsbound2(qp.ncut1());
			D2Leaf<long, std::vector<std::vector<double>>, AngleCrd> d2leaf(
					refls_.get(), 1000000, bnd);
			lstree_.gettree().findneighbor(ls1, d2leaf, bnd);
			std::vector<std::pair<long, double>> &neighbors =
					d2leaf.nnearest().neighbors();
			s1 = 0.0;
			for (auto &n : neighbors) {
				if (n.second < 0.000001)
					continue;  //ignore self
				const std::vector<double> & temp = refls_->at(n.first);
				s1 += matchls(ls1, temp, qp.ncut1());
			}
			s1 /= refls_->size();
		} while (s1 < 0.15 && qp.ncut1() < 6);
	} else {
		double bnd = lsbound2(qp.ncut1());
		D2Leaf<long, std::vector<std::vector<double>>, AngleCrd> d2leaf(
				refls_.get(), 1000000, bnd);
		lstree_.gettree().findneighbor(ls1, d2leaf, bnd);
		std::vector<std::pair<long, double>> &neighbors =
				d2leaf.nnearest().neighbors();
		s1 = 0.0;
		for (auto &n : neighbors) {
			if (n.second < 0.000001)
				continue;  //ignore self
			const std::vector<double> & temp = refls_->at(n.first);
			s1 += matchls(ls1, temp, qp.ncut1());
		}
		s1 /= refls_->size();
	}
	double s2;
	if (qp.ncut2() < 0) {
		qp.ncut2() = -1;
		do {
			++(qp.ncut2());
			double bnd = lsbound2(qp.ncut2());
			D2Leaf<long, std::vector<std::vector<double>>, AngleCrd> d2leaf(
					refls_.get(), 1000000, bnd);
			lstree_.gettree().findneighbor(ls2, d2leaf, bnd);
			std::vector<std::pair<long, double>> &neighbors =
					d2leaf.nnearest().neighbors();
			s2 = 0.0;
			for (auto &n : neighbors) {
				if (n.second < 0.000001)
					continue;  //ignore self
				const std::vector<double> & temp = refls_->at(n.first);
				s2 += matchls(ls2, temp, qp.ncut2());
			}
			s2 /= refls_->size();
		} while (s2 < 0.15 && qp.ncut2() < 6);
	} else {
		double bnd = lsbound2(qp.ncut2());
		D2Leaf<long, std::vector<std::vector<double>>, AngleCrd> d2leaf(
				refls_.get(), 1000000, bnd);
		lstree_.gettree().findneighbor(ls2, d2leaf, bnd);
		std::vector<std::pair<long, double>> &neighbors =
				d2leaf.nnearest().neighbors();
		s2 = 0.0;
		for (auto &n : neighbors) {
			if (n.second < 0.000001)
				continue;  //ignore self
			const std::vector<double> & temp = refls_->at(n.first);
			s2 += matchls(ls2, temp, qp.ncut2());
		}
		s2 /= refls_->size();
	}
	return s1 * s2;
}
void ScoreSitePair::score(QueryPair &qp,bool silent) {
	QueryPair swapped = qp.swapped();
	double srx;
	if (qp.alpha() < 0.0) {
		qp.alpha() = ALPHA0 / 1.2;
		int ntry = 0;
		do {
			qp.alpha() *= 1.2;
			swapped.alpha() = qp.alpha();
			srx = 0.5 * (refdmneighborsum(qp) + refdmneighborsum(swapped));
		} while (srx < 5.0 && ++ntry < 4);
	} else {
		srx = 0.5 * (refdmneighborsum(qp) + refdmneighborsum(swapped));
	}
	double srls = reflsneighborsum(qp);
	swapped.ncut2()=qp.ncut1();
	swapped.ncut1()=qp.ncut2();
//	std::cout <<"---- " <<*srx<<" ---- "<<srls<<"  "<<nlscut1<<"  "<<nlscut2<<std::endl;
//	*srx*=srls;
	double sx = 0.5 * (neighborsum(qp) + neighborsum(swapped));
	std::vector<double> sc(2);
	sc[0] = sx;
	sc[1] = srx * srls;
	qp.scores().push_back(sc);
	if(!silent) std::cout <<sc[0]<<" "<<sc[1]<<std::endl;
}
std::vector<long> NSPpdbstatistics::splitsites(
		const std::vector<BackBoneSite> &sites, int ngroups) {
	long grpsize = sites.size() / ngroups;
	assert(grpsize > 50);
	std::vector<long> grpbegins;
	grpbegins.push_back(0);
	for (int i = 1; i < ngroups; ++i) {
		int guess = i * grpsize - 50;
		std::string pdbidold = (sites.begin() + guess)->pdbid;
		while (guess < sites.size()
				&& (sites.begin() + guess)->pdbid == pdbidold)
			++guess;
		grpbegins.push_back(guess);
	}
	grpbegins.push_back(sites.size());
	return grpbegins;
}
void NSPpdbstatistics::drawquerypairs(
		const std::vector<NSPproteinrep::BackBoneSite> &sites, double dcamin,
		double dcamax, int nskip, int nsample, int nrandomtime,
		std::vector<QueryPair> *qp_native, std::vector<QueryPair> *qp_random) {
	int minsep = 5;
	auto start2 = sites.begin() + 1;
	long npairs = 0;
	for (auto iter1 = start2 + minsep; iter1 != sites.end() - 1; ++iter1) {
		if (chainstartsite(iter1)) {
//			std::cout<<"Number of chains: " <<++nchain <<" Number of Pairs: " <<paircount <<std::endl;
			start2 = iter1 + 1;
		}
		if (iter1 - start2 < minsep)
			continue;
		if (!fragstartsite(iter1 - 1, sites.end(), 3))
			continue;
		int sslength = 0;
		if (iter1->sscodechar() != 'C') {
			while (true) {
				if (iter1 - sslength == start2)
					break;
				if ((iter1 - (sslength + 1))->sscodechar()
						== iter1->sscodechar()) {
					++sslength;
					continue;
				} else {
					break;
				}
			}
		}
		BackBoneSite s1 = *iter1;
		int iter2end = sslength > minsep ? sslength : minsep;
		for (auto iter2 = start2; iter2 != iter1 - iter2end; ++iter2) {
			if (!fragstartsite(iter2 - 1, sites.end(), 3))
				continue;
			BackBoneSite s2 = *iter2;
			double dca2 = (s2.cacrd() - s1.cacrd()).squarednorm();
			if (dca2 < dcamin * dcamin || dca2 >= dcamax * dcamax)
				continue;
			++npairs;
			if (npairs <= nskip)
				continue;
			SiteItPair sip = std::make_pair(iter1, iter2);
			qp_native->push_back(QueryPair(sip));
			std::vector<BackBoneSite> rpair = makerandompairconf(sip,
					nrandomtime, dcamin, dcamax);
			for (int n = 0; n < nrandomtime; ++n) {
				SiteItPair rsip = std::make_pair(rpair.cbegin() + 1,
						rpair.cbegin() + 3 * n + 4);
				qp_random->push_back(QueryPair(rsip));
			}
			if (npairs >= nskip + nsample)
				return;
		}
	}
}


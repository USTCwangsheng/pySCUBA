/*
 * estimatecbpairene.cpp
 *
 *  Created on: 2018年1月24日
 *      Author: hyliu
 */

#include "backbone/backbonesite.h"
#include "dstl/randomengine.h"
#include "geometry/rotation.h"

using namespace NSPproteinrep;
struct DisDistr {
	double step;
	int nsteps;
	std::vector<double> densities;
	DisDistr(double istep = 0, int insteps = 0) :
			step(istep), nsteps(insteps) {
		;
	}
	void estimate(const std::vector<double> & rsamples,
			const std::vector<double> &weights) {
		std::vector<double> count(nsteps + 1, 0.0);
		double count_tot = 0.0;
		int sidx=-1;
		for (auto r : rsamples) {
			sidx++;
			int idx = index(r);
			if (idx < nsteps) {
				count[idx] += weights[sidx];
				count_tot += weights[sidx];
			}
		}
		double count_even = count_tot / (double) nsteps;
		densities.resize(nsteps+1,0.0);
		for (int i = 0; i < nsteps; i++) {
			densities[i] = count[i] / count_even;
		}
		densities[nsteps] = 1.0;
	}
	int index(double r) {
		if (r < 0)
			return 0;
		int idx = (int) (r / step);
		if (idx > nsteps)
			idx = nsteps;
		return idx;
	}
	int index(double r) const {
		if (r < 0)
			return 0;
		int idx = (int) (r / step);
		if (idx > nsteps)
			idx = nsteps;
		return idx;
	}
	double bincenter(int index) const {
		return ((double) index + 0.5) * step;
	}
	double densityintpl(double r) const {
		int i1 = index(r - 0.5 * step);
		int i2 = index(r + 0.5 * step);
		if (i1 == i2)
			return densities[i1];
		double r1 = bincenter(i1);
		double r2 = bincenter(i2);
		double alpha = (r - r1) / r2 - r1;
		return densities[i1] * (1.-alpha) + densities[i2] *  alpha;
	}
	/*	double eneintpl(double r) const {
	 int i1=index(r-0.5*step);
	 int i2=index(r+0.5*step);
	 if(i1==i2) return -log(densities[i1]);
	 double r1=bincenter(i1);
	 double r2=bincenter(i2);
	 double alpha=(r-r1)/r2-r1;
	 return -log(densities[i1])*alpha-log(densities[i2])*(1.0-alpha);
	 }*/
};

double shiftangle(double ang) {
	if (ang > 180.0)
		ang -= 360.0;
	else if (ang < -180.0)
		ang += 360.0;
	return ang;
}

class PhiPsiRegion {
public:
	enum {
		TURN, ALPHA, BETA, INTAB
	};
	static int getregionid(double phi, double psi) {
		phi = shiftangle(phi);
		psi = shiftangle(psi);
		if (phi > 0)
			return TURN;
		else if (psi > -100.0 && psi < 30.0)
			return ALPHA;
		else if (psi > 80.0)
			return BETA;
		return INTAB;
	}
};

std::map<std::pair<int, int>, DisDistr> disdistrs;
std::map<std::pair<int, int>, DisDistr> disdisref;
double b0 = 1.5;
struct CBPair {
	CBPair(NSPproteinrep::BackBoneSite &bs1, NSPproteinrep::BackBoneSite &bs2) {
		int reg1 = PhiPsiRegion::getregionid(bs1.phi(), bs1.psi());
		int reg2 = PhiPsiRegion::getregionid(bs2.phi(), bs2.psi());
		phipsiregions =
				reg1 < reg2 ?
						std::make_pair(reg1, reg2) : std::make_pair(reg2, reg1);
		rcb = sqrt((bs1.cbcrd(b0) - bs2.cbcrd(b0)).squarednorm());
		if(bs1.resname=="GLY"|| bs2.resname=="GLY") wght=0.0000001;
		if(bs1.resname=="ALA"|| bs2.resname=="ALA") wght=0.01;

	}
	std::pair<int, int> phipsiregions;
	double rcb { 0.0 };
	double wght{1.0};
};

double disstep = 0.1;
double nsteps = 80;
void estimatedistrs(const std::vector<CBPair> &samples,
		std::map<std::pair<int, int>, DisDistr> &distrs) {
	std::map<std::pair<int, int>, std::vector<double> > rsamples;
	std::map<std::pair<int, int>, std::vector<double>> rsamples_random;
	std::map<std::pair<int,int>,std::vector<double>> wghts;
	for (auto &s : samples) {
		auto &regid = s.phipsiregions;
		if (rsamples.find(regid) == rsamples.end()) {
			rsamples.insert(std::make_pair(regid, std::vector<double>()));
			wghts.insert(std::make_pair(regid,std::vector<double>()));
		}
		rsamples.at(regid).push_back(s.rcb);
		wghts.at(regid).push_back(s.wght);
	}
	for (auto & rs : rsamples) {
		if (distrs.find(rs.first) == distrs.end()) {
			distrs.insert(std::make_pair(rs.first, DisDistr(disstep, nsteps)));
		}
		distrs.at(rs.first).estimate(rs.second,wghts.at(rs.first));
	}
}

CBPair randomcbpair(NSPproteinrep::BackBoneSite bs1,
		NSPproteinrep::BackBoneSite bs2) {
	bs1.translate(-1.0 * bs1.cacrd());
	NSPgeometry::XYZ trans;
	auto &reng = NSPdstl::RandomEngine<>::getinstance();
	int nget = 0;
	double dcamax = 11.0;
	while (true) {
		trans = NSPgeometry::XYZ(reng.realrng(0.0, 1.0), dcamax);
		trans = trans - bs2.cacrd();
		bs2.translate(trans);
		NSPgeometry::XYZ axis(reng.realrng(), 1.0);
		double angle = reng.realrng(0.0, 1.0)() * 180.0;
		NSPgeometry::Rotation rot;
		rot.init(NSPgeometry::QuaternionCrd(axis, angle), bs2.cacrd());
		bs2.rotate(rot);

		bool clashed = atomsclashed(bs1, bs2);
		if (!clashed)
			break;
	}
	bs1.resname="NNN";
	bs2.resname="NNN";
	return CBPair(bs1, bs2);
}
int maxpairs=100000000;
void collectsamplepairs(const std::vector<BackBoneSite> &sites,
		std::vector<CBPair> &trainingpairs, std::vector<CBPair> &randompairs) {
	int minsep = 5;
	auto start2 = sites.begin() + 1;
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
			if (dca2 > 121.0)
				continue;
	//		double rcb2=(s1.cbcrd()-s2.cbcrd()).squarednorm();
	//		if(rcb2<9.0) std::cout <<sqrt(rcb2)<<std::endl<< s1.toString() <<s2.toString();
			trainingpairs.push_back(CBPair(s1, s2));
			randompairs.push_back(randomcbpair(s1, s2));
		}
		if(trainingpairs.size()>maxpairs) break;
	}
	std::cout << "Number of site pairs: " << trainingpairs.size() << std::endl;
}

int main(int argc, char **argv) {
	std::vector<BackBoneSite> sites;
	readbackbonesites(std::string(argv[1]), sites);
	std::vector<CBPair> trainingpairs;
	std::vector<CBPair> randompairs;
	collectsamplepairs(sites, trainingpairs, randompairs);
	DisDistr totaldistr(disstep,nsteps);
	DisDistr refdistr(disstep,nsteps);
	std::vector<double> rsamples;
	std::vector<double> rrefsamples;
	std::vector<double> wghts;
	std::vector<double> refweights;
	for(auto &cbp: trainingpairs) {
		rsamples.push_back(cbp.rcb);
		wghts.push_back(cbp.wght);
	}
	for(auto &cbp: randompairs) {
		rrefsamples.push_back(cbp.rcb);
		refweights.push_back(cbp.wght);
	}
	totaldistr.estimate(rsamples,wghts);
	refdistr.estimate(rrefsamples,refweights);
	std::string filename = "cbrdist.dat";
	std::ofstream ofs;
	ofs.open(filename.c_str());
	for (int i = 0; i <= nsteps; ++i) {
		double r = totaldistr.bincenter(i);
		double di=totaldistr.densities[i];
		if(di>0.000001){
			ofs << r << " " << -log(di/ refdistr.densities[i])
				<< " " << di << " "
				<< refdistr.densities[i] << std::endl;
		} else{
			ofs<<r <<" "<<20.0
					<< " " << di << " "
										<< refdistr.densities[i] << std::endl;
		}
	}
	ofs.close();
/*	estimatedistrs(trainingpairs, disdistrs);
	estimatedistrs(randompairs, disdisref);
	std::vector<std::string> regionnames { "turn", "alpha", "beta", "intab" };
	for (int i = 0; i < 4; i++) {
		for (int j = i; j < 4; j++) {
			std::string filename = regionnames[i] + regionnames[j]
					+ "cbrdist.dat";
			std::ofstream ofs;
			ofs.open(filename.c_str());
			DisDistr & distr = disdistrs.at(std::make_pair(i, j));
			DisDistr &refdistr = disdisref.at(std::make_pair(i, j));
			for (int i = 0; i <= nsteps; ++i) {
				double r = distr.bincenter(i);
				double di=distr.densities[i];
				if(di>0.000001){
					ofs << r << " " << di/ refdistr.densities[i]
						<< " " << di << " "
						<< refdistr.densities[i] << std::endl;
				} else{
					ofs<<r <<" "<<0.0
							<< " " << di << " "
												<< refdistr.densities[i] << std::endl;
				}
			}
			ofs.close();
		}
	}*/
}


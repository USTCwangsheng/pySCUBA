/*
 * testvirtualsidechain.cpp
 *
 *  Created on: 2018年2月3日
 *      Author: hyliu
 */

#include "pdbstatistics/virtualsidechain.h"
#include "dstl/randomengine.h"

using namespace NSPproteinrep;
using namespace NSPpdbstatistics;
using namespace NSPgeometry;
std::vector<std::string> mcatoms{"N","CA","C","O"};
std::vector<std::string> vscatoms{"VB","VG1","VG2","VD1","VD2"};
std::map<std::string,std::vector<double>> distances{
	{"NVB",std::vector<double>()},
	{"CAVB",std::vector<double>()},
	{"CVB",std::vector<double>()},
	{"OVB",std::vector<double>()},
	{"NVG1",std::vector<double>()},
	{"CAVG1",std::vector<double>()},
	{"CVG1",std::vector<double>()},
	{"OVG1",std::vector<double>()},
	{"NVG2",std::vector<double>()},
	{"CAVG2",std::vector<double>()},
	{"CVG2",std::vector<double>()},
	{"OVG2",std::vector<double>()},
	{"NVD1",std::vector<double>()},
	{"CAVD1",std::vector<double>()},
	{"CVD1",std::vector<double>()},
	{"OVD1",std::vector<double>()},
	{"NVD2",std::vector<double>()},
	{"CAVD2",std::vector<double>()},
	{"CVD2",std::vector<double>()},
	{"OVD2",std::vector<double>()},
	{"VBVB",std::vector<double>()},
	{"VBVG1",std::vector<double>()},
	{"VBVG2",std::vector<double>()},
	{"VBVD1",std::vector<double>()},
	{"VBVD2",std::vector<double>()},
	{"VG1VG1",std::vector<double>()},
	{"VG1VG2",std::vector<double>()},
	{"VG1VD1",std::vector<double>()},
	{"VG1VD2",std::vector<double>()},
	{"VG2VG2",std::vector<double>()},
	{"VG2VD1",std::vector<double>()},
	{"VG2VD2",std::vector<double>()},
	{"VD1VD1",std::vector<double>()},
	{"VD1VD2",std::vector<double>()},
	{"VD2VD2",std::vector<double>()}
};
std::map<std::string,std::vector<double>> refdistances{
	{"NVB",std::vector<double>()},
	{"CAVB",std::vector<double>()},
	{"CVB",std::vector<double>()},
	{"OVB",std::vector<double>()},
	{"NVG1",std::vector<double>()},
	{"CAVG1",std::vector<double>()},
	{"CVG1",std::vector<double>()},
	{"OVG1",std::vector<double>()},
	{"NVG2",std::vector<double>()},
	{"CAVG2",std::vector<double>()},
	{"CVG2",std::vector<double>()},
	{"OVG2",std::vector<double>()},
	{"NVD1",std::vector<double>()},
	{"CAVD1",std::vector<double>()},
	{"CVD1",std::vector<double>()},
	{"OVD1",std::vector<double>()},
	{"NVD2",std::vector<double>()},
	{"CAVD2",std::vector<double>()},
	{"CVD2",std::vector<double>()},
	{"OVD2",std::vector<double>()},
	{"VBVB",std::vector<double>()},
	{"VBVG1",std::vector<double>()},
	{"VBVG2",std::vector<double>()},
	{"VBVD1",std::vector<double>()},
	{"VBVD2",std::vector<double>()},
	{"VG1VG1",std::vector<double>()},
	{"VG1VG2",std::vector<double>()},
	{"VG1VD1",std::vector<double>()},
	{"VG1VD2",std::vector<double>()},
	{"VG2VG2",std::vector<double>()},
	{"VG2VD1",std::vector<double>()},
	{"VG2VD2",std::vector<double>()},
	{"VD1VD1",std::vector<double>()},
	{"VD1VD2",std::vector<double>()},
	{"VD2VD2",std::vector<double>()}
};

BackBoneSite  randommoved(NSPproteinrep::BackBoneSite bs1,
		NSPproteinrep::BackBoneSite bs2) {
	NSPgeometry::XYZ trans;
	auto &reng = NSPdstl::RandomEngine<>::getinstance();
	int nget = 0;
	double dcamax = 13.0;
	while (true) {
		trans = NSPgeometry::XYZ(reng.realrng(0.0, 1.0), dcamax);
		trans = trans - bs2.cacrd()+bs1.cacrd();
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
	return bs2;
}
void adddistances(const BackBoneSite &bs1, const VirtualSideChain & vs1,
		const BackBoneSite &bs2, const VirtualSideChain &vs2,bool doref=false){
		double rca2cut=169;
		double rca2=(bs1.cacrd()-bs2.cacrd()).squarednorm();
		if(rca2>=rca2cut) return;
		std::vector<XYZ> bcrd1,bcrd2;
		bs1.getcrd(bcrd1);
		bs2.getcrd(bcrd2);
		std::map<std::string,std::vector<double>> *diss=&distances;
		if(doref) diss=&refdistances;
		for(int im=0;im<4;++im){
			for(int is=0;is<5;++is){
				std::string entry=mcatoms[im]+vscatoms[is];
				double d=distance(bcrd1[im],vs2.getcrd(vscatoms[is]));
				(*diss)[entry].push_back(d);
				d=distance(bcrd2[im],vs1.getcrd(vscatoms[is]));
				(*diss)[entry].push_back(d);
			}
		}
		for(int is1=0;is1<5;++is1){
			for(int is2=0;is2<5;++is2){
				std::string entry=vscatoms[is1]+vscatoms[is2];
				if(diss->find(entry) == diss->end()) entry=vscatoms[is2]+vscatoms[is1];
				double d=distance(vs1.getcrd(vscatoms[is1]),vs2.getcrd(vscatoms[is2]));
				(*diss)[entry].push_back(d);
				d=distance(vs1.getcrd(vscatoms[is2]),vs2.getcrd(vscatoms[is1]));
				(*diss)[entry].push_back(d);
			}
		}
}
struct DisDistr {
	double step;
	int nsteps;
	std::vector<double> densities;
	DisDistr(double istep = 0, int insteps = 0) :
			step(istep), nsteps(insteps) {
		;
	}
	void estimate(const std::vector<double> & rsamples,
			std::vector<double> weights=std::vector<double>()) {
		std::vector<double> count(nsteps + 1, 0.0);
		double count_tot = 0.0;
		int sidx=-1;
		if(weights.empty()){
			weights.assign(rsamples.size(),1.0);
		}
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
};

int main(int argc,char **argv){
	std::vector<BackBoneSite>sites;
	readbackbonesites(std::string(argv[1]), sites);
	int minsep = 5;
	auto start2 = sites.begin();
	long npairs = 0;
	for (auto iter1 = start2 + minsep; iter1 != sites.end() - 1; ++iter1) {
			if (chainstartsite(iter1)) {
	//			std::cout<<"Number of chains: " <<++nchain <<" Number of Pairs: " <<paircount <<std::endl;
				start2 = iter1 + 1;
			}
			if (iter1 - start2 < minsep)
				continue;
			BackBoneSite s1 = *iter1;
			int iter2end = minsep;
			VirtualSideChain vs1(s1.cacrd(),s1.ncrd(),s1.ccrd());
			for (auto iter2 = start2; iter2 != iter1 - iter2end; ++iter2) {
				BackBoneSite s2 = *iter2;
				VirtualSideChain vs2(s2.cacrd(),s2.ncrd(),s2.ccrd());
				adddistances(s1,vs1,s2,vs2);
				BackBoneSite rs2=randommoved(s1,s2);
				VirtualSideChain rvs2(rs2.cacrd(),rs2.ncrd(),rs2.ccrd());
				adddistances(s1,vs1,rs2,rvs2,true);
				++npairs;
				if(npairs %10000 ==0 ){
					std::cout<<npairs<<std::endl;
				}
			}
			if(npairs>=1800000) break;
	}
	for(auto &d:distances){
		std::string filename=d.first+"_distances.dat";
		std::ofstream ofs;
		ofs.open(filename.c_str());
		DisDistr distr(0.1,180);
		DisDistr refdistr(0.1,180);
		distr.estimate(d.second);
		refdistr.estimate(refdistances.at(d.first));
		for (int i = 0; i <= 180; ++i) {
			double r = distr.bincenter(i);
			double di=distr.densities[i];
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
	}
/*	std::vector<XYZ> localcrd(5,XYZ());
//	double rad=180.0/3.14159265;
	for(auto & s:sites){
		VirtualSideChain vsc(s.cacrd(),s.ncrd(),s.ccrd());
		XYZ rca=s.cacrd();
		XYZ rn=s.ncrd();
		XYZ rc=s.ccrd();
		std::vector<XYZ> dvbdca;
		for(int k=0;k<3;++k){
			rc[k] +=0.00005;
			VirtualSideChain vscp(rca,rn,rc);
			rc[k] -=0.0001;
			VirtualSideChain vscm(rca,rn,rc);
			dvbdca.push_back((vscp.getcrd("VB")-vscm.getcrd("VB"))/0.0001);
			rc[k] +=0.00005;
		}
		for(int m=0;m<3;++m){
			for(int k=0;k<3;++k){
				std::cout << dvbdca[k][m]<<" ";
			}
			XYZ dvdx(0,0,0);
			dvdx[m]=1.0;
			std::cout <<"  ::::   "<<
					vsc.redistributederiv("VB",dvdx)[2].toString()<<std::endl;
		}

		for(int a=0;a<5;++a){
			localcrd[a]=localcrd[a]+vsc.getlocalcrd(vscatoms[a]);
		}
		double rab=distance(vsc.getcrd("VB"),s.cacrd());
		double abac=angle(vsc.getcrd("VB"),s.cacrd(),s.ccrd())*rad;
		double dbacn=torsion(vsc.getcrd("VB"),s.cacrd(),s.ccrd(),s.ncrd())*rad;
		std::cout <<"VB: " << rab<<" "<<abac<< " "<< dbacn <<std::endl;
		rab=distance(vsc.getcrd("VG1"),vsc.getcrd("VB"));
		abac=angle(vsc.getcrd("VG1"),vsc.getcrd("VB"),s.cacrd())*rad;
		dbacn=torsion(vsc.getcrd("VG1"),vsc.getcrd("VB"),s.cacrd(),s.ncrd())*rad;
		std::cout <<"VG1: " << rab<<" "<<abac<< " "<< dbacn <<std::endl;
		rab=distance(vsc.getcrd("VG2"),vsc.getcrd("VB"));
		abac=angle(vsc.getcrd("VG2"),vsc.getcrd("VB"),s.cacrd())*rad;
		dbacn=torsion(vsc.getcrd("VG2"),vsc.getcrd("VB"),s.cacrd(),s.ncrd())*rad;
		std::cout <<"VG2: " << rab<<" "<<abac<< " "<< dbacn <<std::endl;
		rab=distance(vsc.getcrd("VD1"),vsc.getcrd("VG1"));
		abac=angle(vsc.getcrd("VD1"),vsc.getcrd("VG1"),vsc.getcrd("VB"))*rad;
		dbacn=torsion(vsc.getcrd("VD1"),vsc.getcrd("VG1"),vsc.getcrd("VB"),s.cacrd())*rad;
		std::cout <<"VD1: " << rab<<" "<<abac<< " "<< dbacn <<std::endl;
		rab=distance(vsc.getcrd("VD2"),vsc.getcrd("VG2"));
		abac=angle(vsc.getcrd("VD2"),vsc.getcrd("VG2"),vsc.getcrd("VB"))*rad;
		dbacn=torsion(vsc.getcrd("VD2"),vsc.getcrd("VG2"),vsc.getcrd("VB"),s.cacrd())*rad;
		std::cout <<"VD2: " << rab<<" "<<abac<< " "<< dbacn <<std::endl;
	}*/
/*
	std::cout<<"Averages: "<<std::endl;
	for(int a=0;a<5;++a){
		std::cout<<vscatoms[a]<<": "<< (localcrd[a]/(double) sites.size()).toString()<<std::endl;
	}*/
}


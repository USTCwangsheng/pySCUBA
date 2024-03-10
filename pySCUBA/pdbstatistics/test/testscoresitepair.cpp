/*
 * testscoresitepair.cpp
 *
 *  Created on: 2017年12月5日
 *      Author: hyliu
 */

#include "pdbstatistics/scoresitepair.h"
#include "dstl/randomengine.h"
#include <thread>
#include <memory>
using namespace NSPproteinrep;
using namespace NSPpdbstatistics;
int nthread;
int ngroups;
std::vector<long> grpbegins;
int nskip;
int nsample;
int nrandomtime;
std::vector<QueryPair> qp_native;
std::vector<QueryPair> qp_random;
std::vector<std::pair<double,double>> disranges{
	{0,6.1},{4.9,7.1},{5.9,8.1},{6.9,9.1},{7.9,10.1},{8.9,11.1},{9.9,12.1}
};
std::vector<double> disratios{1.57,1.23,0.86,0.83,0.96,1.05,1.02};

struct ThreadArg {
	int nthread;
	int id{-1};
	ScoreSitePair *scoresitepair;
	int disrangeid;
};

void threadrunner(ThreadArg *arg){
	std::cout<<"Started thread "<<arg->id <<std::endl;
	int idx=0;
	int ntest=0;
	ScoreSitePair *scoresp=arg->scoresitepair;
/*
	std::ofstream ofs_n;
	std::ofstream ofs_r;
	std::string nativefilename="nativepair_"+std::to_string(arg->disrangeid)+"_thread_"+
			std::to_string(arg->id)+".dat";
	std::string randomfilename="randompair_"+std::to_string(arg->disrangeid)+"_thread_"+
			std::to_string(arg->id)+".dat";
	ofs_n.open(nativefilename.c_str());
	ofs_r.open(randomfilename.c_str());*/
	bool silent{true};
	if(arg->id==0) silent=false;
	for(auto &qp:qp_native){
//		assert(scoresp->covered(sqrt(qp.dca2())));
		int myid=idx%(arg->nthread);
		if(myid == arg->id) {
			scoresp->score(qp,silent);
		}
			++idx;
	}
	idx=0;
	for(auto &qp:qp_random){
//		assert(scoresp->covered(sqrt(qp.dca2())));
		int myid=idx%(arg->nthread);
		if(myid == arg->id) {
			scoresp->score(qp,silent);
		}
			++idx;
	}
}
void dogroup(std::vector<BackBoneSite> &sites, int rangeid,double *tmpltsize, double *refsize){
	ScoreSitePair scoresp(disranges[rangeid].first,disranges[rangeid].second);
	scoresp.builddmtree(sites);
	scoresp.buildrefdmtree(1);
	double pinclude=50000.0/(double) sites.size();
	scoresp.buildlstree(sites,pinclude);
	*tmpltsize=(double) (scoresp.tmpltsize());
	*refsize=(double) (scoresp.refsize());
	std::vector<ThreadArg> threadargs(nthread);
	int tid=0;
	for(auto & a:threadargs){
		a.id=tid++;
		a.disrangeid=rangeid;
		a.scoresitepair=&scoresp;
		a.nthread=nthread;
	}
	std::vector<std::shared_ptr<std::thread>> threads(nthread,nullptr);
	for(int id=0; id<nthread;++id){
		threads[id]=std::shared_ptr<std::thread>(new std::thread(threadrunner,&(threadargs[id])));
//		threadrunner(&(threadargs[id]));
	}
	for(int id=0;id<nthread;++id){
		threads[id]->join();
	}
}
/*
std::vector<double> distref(const std::vector<SiteItPair> &confs, double dismax,
		std::vector<std::pair<double,double>> disranges){
	int npairs=500000;
	int ntimes=npairs/confs.size();
	std::vector<double> res(disranges.size(),0.0);
	double total=0.0;
	for(auto &sip:confs){
		std::vector<BackBoneSite> pairs=makerandompairconf(sip,ntimes,0,dismax);
		for(int i=0;i<ntimes;++i){
			total +=1.0;
			double rca=sqrt((pairs[1].cacrd()-pairs[3*i+4].cacrd()).squarednorm());
			for(int m=0;m<disranges.size();++m){
				if(rca>=disranges[m].first &&rca<disranges[m].second) res[m] +=1.0;
			}
		}
	}
	for(auto & r:res) r/=total;
	return res;
}*/
int main(int argc, char ** argv){
	std::vector<BackBoneSite> sites;
	readbackbonesites(std::string(argv[1]), sites);
	int rangeid=std::stoi(std::string(argv[2]));
	ngroups=std::stoi(std::string(argv[3]));
	nskip=std::stoi(std::string(argv[4]));
	nsample=std::stoi(std::string(argv[5]));
	nrandomtime=std::stoi(std::string(argv[6]));
	nthread=std::stoi(std::string(argv[7]));
	if(argc > 8){
		int rseed= std::stoi(std::string(argv[8]));
		NSPdstl::RandomEngine<>::getinstance().reseed(rseed);
	}
	grpbegins=splitsites(sites,ngroups);
	auto range=ScoreSitePair::coverrange(disranges[rangeid].first,disranges[rangeid].second);
	drawquerypairs(sites,range.first,range.second,nskip,nsample,nrandomtime,
			&qp_native,&qp_random);
	std::vector<double> tmpltsizes(ngroups);
	std::vector<double> rfsizes(ngroups);
	for(int i=0;i<ngroups;++i){
		std::vector<BackBoneSite> subset(grpbegins[i+1]-grpbegins[i]);
		std::copy(sites.begin()+grpbegins[i],sites.begin()+grpbegins[i+1],subset.begin());
		dogroup(subset,rangeid,&tmpltsizes[i],&rfsizes[i]);
	}
	std::ofstream ofs_n;
	std::string nativefilename="nativepair_distange"+std::to_string(rangeid)+".dat";
	ofs_n.open(nativefilename.c_str());
	for(auto &qp:qp_native){
		std::vector<double> scores=qp.summscore(tmpltsizes,rfsizes);
		ofs_n<<"Energy: " << -log(disratios[rangeid]*scores[0]);
		for(auto s:scores) ofs_n<<" "<<s;
		ofs_n<<std::endl;
		qp.saveconf(ofs_n);
	}
	ofs_n.close();
	std::ofstream ofs_r;
	std::string randomfilename="randompair_distange"+std::to_string(rangeid)+".dat";
	ofs_r.open(randomfilename.c_str());
	for(auto &qp:qp_random){
		std::vector<double> scores=qp.summscore(tmpltsizes,rfsizes);
		ofs_r<<"Energy: " << -log(disratios[rangeid]*scores[0]);
		for(auto s:scores) ofs_r<<" "<<s;
		ofs_r<<std::endl;
		qp.saveconf(ofs_r);
	}
	ofs_r.close();
}


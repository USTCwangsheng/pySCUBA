/*
 * pocketdocker.cpp
 *
 *  Created on: 2020年7月23日
 *      Author: hyiu
 */
#include "sampling/pocketdocker.h"
#include "dstl/randomengine.h"
using namespace NSPgeometry;
using namespace NSPintrct;
using namespace NSPproteinrep;
using namespace NSPsampling;
void NbrFndrGrid::setgrid(const NSPgeometry::XYZ &center,
		const std::array<int,3> &sizes, double step) {
	grid=RectGrid(center,sizes,step);
}
void NbrFndrGrid::setgridneighbors(const NSPintrct::IntrctMol & imol,double rcut2){
	cutoff2=rcut2;
	const std::vector<XYZ> & crds=imol.recollectcrds_all();
	gridneighbors=grid.neighborxyzs(crds,rcut2);
}
double BckBnStMatcher::sitedistance2(const NSPintrct::IntrctBlck &iblck,
			const NSPproteinrep::AAConformer &aaconf){
	static std::array<std::string,3> bckatms{"N","CA","C"};
	std::vector<XYZ> crd1;
	std::vector<XYZ> crd2;
	for(auto &atnm:bckatms) {
		crd1.push_back(iblck.getcrd(atnm));
		crd2.push_back(aaconf.getglobalcrd().at(atnm));
	}
	crd1.push_back(InternaltoXYZ(crd1[1],crd1[2],crd1[0],3.5,1.911,2.094));
	crd2.push_back(InternaltoXYZ(crd2[1],crd2[2],crd2[0],3.5,1.911,2.094));
	double dist2=0.0;
	for(int i=0;i<4;i++) dist2 +=(crd1[i]-crd2[i]).squarednorm();
	return dist2;
}
NSPdstl::Idx2D BckBnStMatcher::matchsite(const NSPintrct::IntrctMol & imol,
		const NSPproteinrep::AAConformer &aaconf,double &score) const{
		NSPdstl::Idx2D res;
		score=-1.0;
		for (auto &blkidx:candidatesites){
			double dist2=sitedistance2(imol.getblck(blkidx),aaconf);
			if(dist2<maxdist2){
				 res=blkidx;
				 double alpha=dist2/maxdist2;
				 score=1-alpha*alpha*alpha;
			}
		}
		return res;
}
void PocketMover::stepback(NSPproteinrep::PocketModeler &pm){
	pm.copycrdfrom(this->crd_before);
	pocketclashes=pm.totalclashes();
}

bool PocketMover::randommv(NSPproteinrep::PocketModeler &pm){
	auto &rng=NSPdstl::RandomEngine<>::getinstance().realrng(0,1);
	movencomponents.clear();
	crd_before=pm.copycrdto();
	if(rng()<prigidmv){
		pm.rigidmv(maxstep_trans,maxstep_rotate,ligandcentersphere);
		for(auto &c:pm.pocket().components()) movencomponents.insert(c.first);
	} else {
		int rsize=pm.joints().size()+pm.flexrotamers().size();
		auto irng=NSPdstl::RandomEngine<>::getinstance().intrng(0,rsize-1);
			int r=irng();
			if(r <pm.joints().size()) {
				pm.mvaroundjoint(pm.joints().at(r),pocketclashes);
				movencomponents=pm.joints().at(r).mvcmpnts;
			} else{
				auto &fr=pm.flexrotamers().at(r-pm.joints().size());
				pm.mvrotamer(fr,pocketclashes);
				movencomponents.insert(fr.cmpntid);
			}
	}
}
bool PocketMover::mvtosphere(NSPproteinrep::PocketModeler &pm){
	auto &rng=NSPdstl::RandomEngine<>::getinstance().realrng(0,1);
	crd_before=pm.copycrdto();
	movencomponents.clear();
	XYZ nc=ligandcentersphere.randompoint(rng);
	XYZ oc=pm.ligandcenter();
	pm.rigidtranslation(nc-oc);
	for(auto &c:pm.pocket().components()) movencomponents.insert(c.first);
	return true;
}
bool PocketMover::mvtofit(NSPproteinrep::PocketModeler &pm){
	assert(refcrd.size()>2);
	std::vector<XYZ> crds,crdsref;
	crd_before=pm.copycrdto();
	for(auto & c:refcrd){
		crdsref.push_back(c.second);
		crds.push_back(pm.pocket().components().at(c.first.cmpntID).getglobalcrd().at(c.first.atmnm));
	}
	NSPgeometry::QuatFit qf;
	qf.setup(crdsref,crds);
	auto rt=qf.getRigidTransform();
	for(auto &c:pm.pocket().components()){
		for (auto &x:c.second.getglobalcrd()) rt.apply(&x.second);
	}
	movencomponents.clear();
	for(auto &c:pm.pocket().components()) movencomponents.insert(c.first);
	return true;
}
double clashscore(const AtomIP &sa1, const AtomIP &sa2, double dist2 ){
	static const double fac = pow(2, 1 / 6.0);
	double sigma;
	if ((sa1.hbdonor && sa2.hbacceptor) || (sa1.hbacceptor && sa2.hbdonor)) {
		sigma = 0.5 * (sa1.sigmahb + sa2.sigmahb);
	} else if ((sa1.hbdonor && sa2.hbdonor)
			|| (sa1.hbacceptor && sa2.hbacceptor)) {
		sigma = 0.55 * (sa1.sigma + sa2.sigma);
	} else {
		sigma = 0.5 * (sa1.sigma + sa2.sigma);
	}
	double rcut2 = sigma * sigma * 1.2;
	double sl = fac * sigma;
	double eclash;
	if (dist2 > sl * sl) {
		eclash=0.0;
	}else {
		double r6 = dist2*dist2*dist2;
		double r12 = r6 * r6;
		double sigma3 = sigma * sigma * sigma;
		double sigma6 = sigma3 * sigma3;
		double sigma12 = sigma6 * sigma6;
		eclash= 4.0 *  (sigma12 / r12 - sigma6 / r6) +1.0;
	}
	if(eclash>500.0) eclash=500.0;
	return eclash/500.0;
}
double PocketDocker::score_clashes(const NSPproteinrep::IsolatedPocket::ComponentID &comp){
	auto &aaconf=pocketmodeler_->pocket().components().at(comp);
	const BlckTopo *aatopo;
	if(aaconf.isnaturalaa()) aatopo=&(BlckTopo::getblcktopo_std(aaconf.residuename));
	else aatopo=&(BlckTopo::getblcktopo_nonstd(aaconf.residuename,&aaconf));
	double score=0.0;
	bool isligand=pocketmodeler_->pocket().ligandcomponents().find(comp)
						!=pocketmodeler_->pocket().ligandcomponents().end();
	for(auto &crd:aaconf.getglobalcrd()){
		auto & atomips=aatopo->atomips[aatopo->atomindex(crd.first)];
		auto  cnbrs=nbrfndrgrid_.candidateneighbors(crd.second);
		for(int i:cnbrs){
			auto &presidue=(imol_->getmolsystm())[imol_->residueofatoms()[i]];
			if(presidue.regioncode !=IntrctBlck::PROTEIN_CORE  & !isligand) continue;
			int iatom=i-presidue.aoffset;
			auto & crd2=presidue.crds[i];
			double dist2=(crd2-crd.second).squarednorm();
			if(dist2 >= nbrfndrgrid_.cutoff2) continue;
			if(isligand &&
					presidue.lexcludedatoms.find(iatom) !=presidue.lexcludedatoms.end()) continue;
			auto &atomips2=presidue.topo->atomips[i];
			score+=clashscore(atomips, atomips2,dist2);
		}
	}
	return score;
}
double PocketDocker::score_sitematch(const NSPproteinrep::IsolatedPocket::ComponentID &comp){
	auto &aaconf=pocketmodeler_->pocket().components().at(comp);
	double score;
	auto blk=bckbnstmatcher_.matchsite(*imol_,aaconf,score);
	sitematchscores_.erase(comp);
	if(score >0){
		sitematchscores_[comp]= std::make_pair(blk,score);
	}
	return score;
}
double PocketDocker::score_clashes(){
	for(auto &comp:pocketmover_.movencomponents){
		clash_scores_[comp]=score_clashes(comp);
	}
	totalclashscore_=0.0;
	for(auto & cs:clash_scores_) totalclashscore_ +=cs.second;
	return totalclashscore_;
}
double PocketDocker::score_sitematch(){
	for(auto &comp:pocketmover_.movencomponents){
			score_sitematch(comp);
		}
		 totalmatchscore_=0.0;
		for(auto & cs:sitematchscores_) totalmatchscore_ +=cs.second.second;
		return totalmatchscore_;
}
bool PocketDocker::micromove(){
	double old_clashscore=totalclashscore_;
	double old_matchscore=totalmatchscore_;
	auto save_c=this->clash_scores_;
	auto save_m=this->sitematchscores_;
	int ntrial=0;
#define MAXTRIAL 1000
	bool moven=false;
	while (ntrial <MAXTRIAL){
		pocketmover_.randommv(*pocketmodeler_);
		score_clashes();
		score_sitematch();
		if(totalclashscore_>=old_clashscore && totalmatchscore_>=old_matchscore){
			break;
			moven=true;
		}
		pocketmover_.stepback(*pocketmodeler_);
		this->clash_scores_=save_c;
		this->sitematchscores_=save_m;
		++ntrial;
	}
	return moven;
}


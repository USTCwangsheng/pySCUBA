/*
 * molmodeler.cpp
 *
 *  Created on: 2019年12月19日
 *      Author: hyiu
 */
#include "iblock/molmodeler.h"
#include "designseq/RotamerLib.h"
#include "geometry/calculators.h"
#include "backbone/backbonebuilder.h"
using namespace NSPintrct;
using namespace NSPdesignseq;
using namespace NSPproteinrep;
LoopReservoir MolModeler::mkloopreservoir(const SiteIdx &leftsite, const SiteIdx &rightsite) {
	assert(leftsite.first==rightsite.first);
	BackBoneSite right,left;
	if(rightsite.second<target_->nresidues(rightsite.first)){
		right=(*target_)(rightsite).getbackbonesite();
		right.data_[BackBoneSite::PHI]=this->getphipsi(rightsite).phi;
		right.data_[BackBoneSite::PSI]=this->getphipsi(rightsite).psi;
	}
	if(leftsite.second>=0){
		left=(*target_)(leftsite).getbackbonesite();
			left.data_[BackBoneSite::PHI]=this->getphipsi(leftsite).phi;
			left.data_[BackBoneSite::PSI]=this->getphipsi(leftsite).psi;
	}
	if(leftsite.second==-1){
		LoopReservoir res;
		res.init_nterminus(right);
	} else if(rightsite.second==target_->nresidues(rightsite.first)){
		LoopReservoir res;
		res.init_cterminus(left);
	}
	return LoopReservoir(left,right);
}
const std::vector<MolModeler::SiteIdx> & MolModeler::getnearbysites(
		const SiteIdx & site) {
	if (siteneighbors_.find(site) != siteneighbors_.end())
		return siteneighbors_.at(site);
	auto &blck = (*target_)(site);
	auto & bsinchains = target_->bsinchains;
	NSPgeometry::XYZ rca = blck.getrepcrd();
	siteneighbors_[site] = std::vector<SiteIdx>();
	for (auto &c : bsinchains) {
		for (auto &r : c) {
			double dist =
					(rca - (*target_)(NSPdstl::Idx2D{r.chainid,r.resid}).getrepcrd()).squarednorm();
			if (dist > this->neighbor_cutoff2_)
				continue;
			if (r.chainid == site.first && r.resid == site.second)
				continue;
			siteneighbors_[site].push_back(SiteIdx{r.chainid, r.resid});
		}
	}
	return siteneighbors_.at(site);
}
NSPdesignseq::Phipsi MolModeler::getphipsi(const SiteIdx &site) const {
	assert(phipsis_.find(site) != phipsis_.end());
	return phipsis_.at(site);
}
NSPdesignseq::Phipsi MolModeler::getphipsi(const SiteIdx &site) {
	if (phipsis_.find(site) != phipsis_.end())
		return phipsis_.at(site);
	double phi = 360;
	double psi = 360;
	NSPgeometry::XYZ xr, xn, xca, xc, xl;
	xn = (*target_)(site).getcrd("N");
	xca = (*target_)(site).getcrd("CA");
	xc = (*target_)(site).getcrd("C");
	if (site.second > 0) {
		xr = (*target_)( { site.first, site.second - 1 }).getcrd("C");
		phi = NSPgeometry::torsion(xr, xn, xca, xc) * 180.0 / 3.14159265;
	}
	if (site.second < target_->nresidues(site.first) - 1) {
		xl = (*target_)( { site.first, site.second + 1 }).getcrd("N");
		psi = NSPgeometry::torsion(xn, xca, xc, xl) * 180.0 / 3.14159265;
	}
	phipsis_[site] = Phipsi(phi, psi);
	return phipsis_.at(site);
}
/**
 * access the backbonesites_ data
 *
 *This non-const overload setup the data when needed
 */
const std::vector<std::vector<NSPproteinrep::BackBoneSite>> &MolModeler::getbackbonesites(){
if(!backbonesites_.empty()) return backbonesites_;
backbonesites_.assign(target_->nchains(),std::vector<NSPproteinrep::BackBoneSite>());
for(int c=0;c<backbonesites_.size();++c){
	for(auto & r:target_->getmolsystm().vec(c)){
		if(r.isaminoacid)backbonesites_[c].push_back(r.getbackbonesite());
	}
}
return backbonesites_;
}
std::array<int,2> MolModeler::nclashes(const IntrctBlck & blk1,
		const IntrctBlck &blk2) const {
	std::array<int,2> nc{0,0};
	const IntrctBlck *b1 = &blk1;
	const IntrctBlck *b2 = &blk2;
	if ((blk1.chainseq() > blk2.chainseq())
			|| (blk1.chainseq() == blk2.chainseq()
					&& blk1.resseq() > blk2.resseq())) {
		b1 = &blk2;
		b2 = &blk1;
	}
	std::array<int, 3> nnb;
	auto nbors = b1->getneighbors(*b2, IntrctPara(), nnb);
	if(nnb[0] !=0){
		for(int i=0;i<nbors[0].size();++i){
			for (int j : nbors[0][i]) {
					double d2 = (b1->getcrd(i) - b2->getcrd(j)).squarednorm();
					if (d2 > this->atomic_clash_nhb_)
						continue;
					if (d2 > this->atomic_clash_hb_) {
						if ((b1->topo->atomips[i].hbdonor
								&& b2->topo->atomips[j].hbacceptor)
								|| (b2->topo->atomips[j].hbdonor
										&& b1->topo->atomips[i].hbacceptor))
							continue;
					}
					nc[0]++;
				}
		}
	}
	if (nnb[1] == 0)
		return nc;
	for (int i = 0; i < nbors[1].size(); ++i) {
		for (int j : nbors[1][i]) {
			double d2 = (b1->getcrd(i) - b2->getcrd(j)).squarednorm();
			if (d2 > this->atomic_clash_nhb_)
				continue;
			if (d2 > this->atomic_clash_hb_) {
				if ((b1->topo->atomips[i].hbdonor
						&& b2->topo->atomips[j].hbacceptor)
						|| (b2->topo->atomips[j].hbdonor
								&& b1->topo->atomips[i].hbacceptor))
					continue;
			}
			nc[1]++;
		}
	}
	return nc;
}
const NSPgeometry::LocalFrame & MolModeler::getlocalframe(
		const SiteIdx &site) const {
	assert(sitelocalframes_.find(site) != sitelocalframes_.end());
	return sitelocalframes_.at(site);
}
const NSPgeometry::LocalFrame & MolModeler::getlocalframe(const SiteIdx &site) {
	if (sitelocalframes_.find(site) != sitelocalframes_.end())
		return sitelocalframes_.at(site);
	BackBoneSite bs;
	std::vector<NSPgeometry::XYZ> bcrd;
	const IntrctBlck & blk = (*target_)(site);
	bcrd.push_back(10.0*blk.getcrd("N"));
	bcrd.push_back(10.0*blk.getcrd("CA"));
	bcrd.push_back(10.0*blk.getcrd("C"));
	bcrd.push_back(10.0*blk.getcrd("O"));
	bs.changecrd(bcrd);
	sitelocalframes_[site] = NSPdesignseq::getBackboneSiteLocalFrame(bs);
	return sitelocalframes_.at(site);
}

IntrctBlck MolModeler::mkfloatingblck(const SiteIdx &site, Rotamer *rt) {
	IntrctBlck blk = (*target_)(site);
	int natoms_old = blk.natoms();
	blk.topo = &(BlckTopo::getblcktopo_std(rt->triName));
	if (natoms_old != blk.natoms()) {
		blk.crds.resize(blk.natoms());
		blk.crds[blk.natoms() - 1] = (*target_)(site).getcrd("O");
		blk.crds[blk.natoms() - 2] = (*target_)(site).getcrd("C");
	}
	auto cs = getlocalframe(site);
	std::vector<NSPgeometry::XYZ> sccrds;
	rt->buildSidechain(cs, sccrds);
	for(auto &x:sccrds) x=x*0.1;
	blk.changesccrds(sccrds);
	return blk;
}
IntrctBlck MolModeler::mkfloatingblck(const SiteIdx &site,
		const std::string restype, const NSPproteinrep::BackBoneSite &bs) const {
	IntrctBlck blk;
	blk.topo = &(BlckTopo::getblcktopo_std(restype));
	blk.crds.resize(blk.natoms());
	blk.crds[0]=bs.ncrd()*A2NM;
	blk.crds[1]=bs.cacrd()*A2NM;
	if(restype !="GLY") blk.crds[2]=bs.cbcrd()*A2NM;
	blk.crds[blk.natoms() - 1] = bs.ocrd()*A2NM;
	blk.crds[blk.natoms() - 2] = bs.ccrd()*A2NM;
	blk.setimol(target_,site.first,site.second);
	return blk;
}
std::vector<IntrctBlck> MolModeler::newloop(LoopConfMaker &confmaker,
		const std::vector<std::string> &loopresiduetypes,
		int acceptedclashes, int maxtestnumber, bool helixinloop){
	std::vector<IntrctBlck> savedconf;
	int length=loopresiduetypes.size();
	int ntest=0;
	int nminclash=10000000;
	while(ntest++<maxtestnumber){
		//std::cout << "try " << ntest << "times newloop" << std::endl;
		if (ntest == 1 || ntest == 10 || ntest  == 50 || ntest % 100 == 0) { std::cout << "try " << ntest << "times newloop" << std::endl; }
		//if(ntest%20==0){ std::cout << "try " << ntest << "times newloop"<< std::endl; }
		auto & conf=*(confmaker.reservoir.popconf(length, helixinloop)); // original is: popconf(length)
		//std::cout  << ntest << "times done" << std::endl;
		if (conf.empty()&& ntest== maxtestnumber-1) return savedconf;
		if(conf.empty()) continue;
		//std::cout << "conf.size() cheack"<< conf.size() << std::endl;
		if (conf.size() != length ) continue;
		std::vector<IntrctBlck> blcks;
		int idx=0;
		for(auto &bs:conf){
			blcks.push_back(this->mkfloatingblck(
					SiteIdx{-10,-10},
					loopresiduetypes[idx],bs));
		}
		int nc_mc=0;
		for(auto &blk1:blcks){
			for(int c=0;c<target_->nchains();++c){
				for(int r=0;r<target_->nresidues(c);++r){
					if(c==confmaker.leftsite.first){
						if(r>=confmaker.leftsite.second && r<=confmaker.rightsite.second) continue;
					}
					auto nc=nclashes(blk1,(*target_)({c,r}));
					nc_mc +=nc[0];
				}
			}
		}
		if(nc_mc<=acceptedclashes){
			//std::cout << "nc_mc<=acceptedclashes" << std::endl;
			savedconf=blcks;
			break;
		}
		if(nminclash < nc_mc){
			//std::cout << "nminclash < nc_mc" << std::endl;
			savedconf=blcks;
			nminclash=nc_mc;
		}
	}
	//if (conf.empty()) return savedconf;
	if(!savedconf.empty()){
		//std::cout << "!savedconf.empty()" << std::endl;
		for(int i=0;i<savedconf.size();++i){
			savedconf[i].setimol(target_,confmaker.leftsite.first,confmaker.leftsite.second+i+1);
		}
	}
	return savedconf;
}
//std::vector<IntrctBlck> MolModeler::newloop(LoopConfMaker &confmaker,
//		const std::vector<std::string> &loopresiduetypes,
//		int acceptedclashes, int maxtestnumber){
//	std::vector<IntrctBlck> savedconf;
//	int length=loopresiduetypes.size();
//	int ntest=0;
//	int nminclash=10000000;
//	while(ntest++<maxtestnumber){
//		auto & conf=*(confmaker.reservoir.popconf(length));
//		if(conf.empty()) return savedconf;
//		std::vector<IntrctBlck> blcks;
//		int idx=0;
//		for(auto &bs:conf){
//			blcks.push_back(this->mkfloatingblck(
//					SiteIdx{-10,-10},
//					loopresiduetypes[idx],bs));
//		}
//		int nc_mc=0;
//		for(auto &blk1:blcks){
//			for(int c=0;c<target_->nchains();++c){
//				for(int r=0;r<target_->nresidues(c);++r){
//					if(c==confmaker.leftsite.first){
//						if(r>=confmaker.leftsite.second && r<=confmaker.rightsite.second) continue;
//					}
//					auto nc=nclashes(blk1,(*target_)({c,r}));
//					nc_mc +=nc[0];
//				}
//			}
//		}
//		if(nc_mc<=acceptedclashes){
//			savedconf=blcks;
//			break;
//		}
//		if(nminclash < nc_mc){
//			savedconf=blcks;
//			nminclash=nc_mc;
//		}
//	}
//	if(!savedconf.empty()){
//		for(int i=0;i<savedconf.size();++i){
//			savedconf[i].setimol(target_,confmaker.leftsite.first,confmaker.leftsite.second+i+1);
//		}
//	}
//	return savedconf;
//}
IntrctMol MolModeler::mkimol_newloop
	 	 (const SiteIdx &siteleft, const SiteIdx &siteright,
	 			 const std::vector<IntrctBlck> &newloop) const {
		IntrctMol imol;
		for(int c=0;c<target_->nchains();++c){
			if(siteleft.first!=c){
				imol.addchain(target_->getmolsystm().vec(c));
			} else {
				std::vector<IntrctBlck> nc;
				for(int i=0;i<=siteleft.second;++i){
					nc.push_back((*target_)({c,i}));
				}
				for(auto &r:newloop) nc.push_back(r);
				for(int i=siteright.second;i<target_->nresidues(c);++i){
					nc.push_back((*target_)({c,i}));
				}
				imol.addchain(nc);
			}
		}
		imol.natoms(); //<*call this to make sure atom offsets of intrctblcks in imol_ properly set
		return imol;
	}
IntrctMol MolModeler::mergechains(int firstchain, int secondchain,
		SiteIdx &leftsite, SiteIdx &rightsite) const{
		 IntrctMol imol;
		 std::vector<IntrctBlck> mergedchain=target_->getmolsystm().vec(firstchain);
		 for(auto &b:target_->getmolsystm().vec(secondchain)) mergedchain.push_back(b);
		 int newchainid=(firstchain<secondchain)?firstchain:secondchain;
		 leftsite=SiteIdx(newchainid,target_->nresidues(firstchain)-2);
		 rightsite=SiteIdx(newchainid,target_->nresidues(firstchain)+1);
		 for(int i=0;i<target_->nchains();++i){
			 if(i!=firstchain && i!=secondchain)
				 imol.addchain(target_->getmolsystm().vec(i));
			 else if(i==newchainid) imol.addchain(mergedchain);
		 }
		 imol.natoms(); //this call setup the blk atom index offsets
		 return imol;
	 }
void NSPintrct::ssregions(const std::vector<char> &ssseq,std::vector<std::pair<int,int>> &helixregions,
	std::vector<std::pair<int,int>> &strandregions){
	char ssprev=' ';
	int idx=0;
	int hstart=-1;
	int sstart=-1;
	for(auto ss:ssseq){
		if(ss !=ssprev) {
			if(ssprev=='H'){
				helixregions.push_back(std::make_pair(hstart,idx-hstart));
				hstart=-1;
			} else if(ssprev=='E'){
				strandregions.push_back(std::make_pair(sstart,idx-sstart));
				sstart=-1;
			}
			if(ss=='H') {
				hstart=idx;
				sstart=-1;
			} else if(ss=='E'){
				sstart=idx;
				hstart=-1;
			}
		}
		ssprev=ss;
		idx++;
	}
	if(hstart !=-1){
		helixregions.push_back(std::make_pair(hstart,idx-hstart));
	}
	if(sstart !=-1){
			strandregions.push_back(std::make_pair(sstart,idx-sstart));
		}
}
void MolModeler::buildnewchain( NSPgeometry::XYZ r0, NSPgeometry::XYZ direction,
		const std::vector<char> &ssseq,
		const std::vector<std::string> & aaseq){
	int len=ssseq.size();
	BackBoneSite bs0;
	genbackbonesite(nullptr,false,-120.0,120.0,&bs0);
	 std::vector<std::pair<int,int>>  helixregions;
	std::vector<std::pair<int,int>>  strandregions;
	ssregions(ssseq,helixregions,strandregions);
	std::vector<NSPproteinrep::BackBoneSite> bss=
			BackBoneBuilder::buildforwardbackbone(len,
			bs0,  helixregions,
			strandregions, std::set<int>());
	BackBoneBuilder::movechainto(r0,direction,true,bss);
	std::vector<IntrctBlck> blcks;
	for(auto &bs:bss){
		std::string restype="GLY";
		blcks.push_back(this->mkfloatingblck(
				std::make_pair(-10,-10),
				restype,bs));
	}
	target_->addchain(blcks);
	if(!aaseq.empty()){
		int cid=target_->nchains()-1;
		for(int i=0;i<len;++i){
			this->buildsidechain(std::make_pair(cid,i),aaseq[i],false);
		}
	}
}

// add buildssele by miaoyy
void MolModeler::buildssele(std::vector<char> sseq, std::vector<std::vector<double>> sscrd,std::vector<int> sslen,int loop_option,std::vector<char> &finalss,std::vector<std::string> &ss_info,  std::vector<std::string>& ssEH, std::vector<std::string>& ssE, std::vector<std::string>& ssH, std::vector<std::string>& ssC,std::vector<char> direction){
	std::vector<std::vector<BackBoneSite>> ssele;
	for(int i = 0; i < sseq.size(); i=i+2){
	assert(sseq[i]==sseq[i+1]);
	char sstype=sseq[i];
	XYZ ss_start=XYZ(sscrd[i][0],sscrd[i][1],sscrd[i][2]);
	XYZ ss_end=XYZ(sscrd[i+1][0],sscrd[i+1][1],sscrd[i+1][2]);
	XYZ direc_N=ss_end;
        XYZ direc_C=XYZ(-sscrd[i+1][0],-sscrd[i+1][1],-sscrd[i+1][2]);
  // XYZ direc_N=XYZ(sscrd[i+1][0]-sscrd[i][0],sscrd[i+1][1]-sscrd[i][1],sscrd[i+1][2]-sscrd[i][2]);
	//XYZ direc_C=XYZ(sscrd[i][0]-sscrd[i+1][0],sscrd[i][1]-sscrd[i+1][1],sscrd[i][2]-sscrd[i+1][2]);

	double dis=std::sqrt((ss_start.x_-ss_end.x_)*(ss_start.x_-ss_end.x_)+(ss_start.y_-ss_end.y_)*(ss_start.y_-ss_end.y_)+(ss_start.z_-ss_end.z_)*(ss_start.z_-ss_end.z_));
	int sslength =sslen[int(i/2)];
	int ss_direction=direction[int(i/2)];
  BackBoneSite nsite;
	genbackbonesite(nullptr, false, 0.0,0.0, &nsite);
	if(sstype=='E'){
		std::vector<std::pair<int,int>> strand;
		strand.push_back(std::make_pair(0,sslength));
		if (ss_direction=='N'){
			std::vector<BackBoneSite> chain=BackBoneBuilder::buildforwardbackbone(sslength,nsite, std::vector<std::pair<int,int>>(),strand,std::set<int> ());
			BackBoneBuilder::movechainto(ss_start,direc_N,true,chain);
			ssele.push_back(chain);
		}
		else if (ss_direction == 'C'){
			std::vector<BackBoneSite> chain=BackBoneBuilder::buildbackwardbackbone(sslength,nsite, std::vector<std::pair<int,int>>(),strand,std::set<int> ());
			BackBoneBuilder::movechainto(ss_start,direc_C,false,chain);
			ssele.push_back(chain);
		}

		}
  else if(sstype=='H'){
		std::vector<std::pair<int,int>> helix;
		helix.push_back(std::make_pair(0,sslength));
		if (ss_direction=='N'){
			std::vector<BackBoneSite> chain=BackBoneBuilder::buildforwardbackbone(sslength,nsite, helix ,std::vector<std::pair<int,int>>(),std::set<int> ());
			BackBoneBuilder::movechainto(ss_start,direc_N,true,chain);
			ssele.push_back(chain);
		}
		else if (ss_direction == 'C'){
			std::vector<BackBoneSite> chain=BackBoneBuilder::buildforwardbackbone(sslength,nsite, helix ,std::vector<std::pair<int,int>>(),std::set<int> ());
			BackBoneBuilder::movechainto(ss_start,direc_C,false,chain);
			ssele.push_back(chain);
		}


		}
	}
	std::vector<int> looplen;
	for(int i=0;i<ssele.size()-1;i++){
		// std::cout<<ssele[i].back().cacrd().x_<<ssele[i].back().cacrd().y_<<ssele[i].back().cacrd().z_<<std::endl;
		double eledis=std::sqrt((ssele[i].back().cacrd().x_-ssele[i+1][0].cacrd().x_)*(ssele[i].back().cacrd().x_-ssele[i+1][0].cacrd().x_)
		+(ssele[i].back().cacrd().y_-ssele[i+1][0].cacrd().y_)*(ssele[i].back().cacrd().y_-ssele[i+1][0].cacrd().y_)
		+(ssele[i].back().cacrd().z_-ssele[i+1][0].cacrd().z_)*(ssele[i].back().cacrd().z_-ssele[i+1][0].cacrd().z_));
		int looplength=(int)(eledis/3.75)+2;//+3在直线距离外多加两个，可变最好
		if (looplength<3){looplength=3;}
                looplen.push_back(looplength);

	}

	//link i i+1
	std::vector<std::vector<std::shared_ptr<std::vector<BackBoneSite>>>> loops;
	if (loop_option)
	{	for(int i=0;i<ssele.size()-1;i++){
			int maxtrytime=0;
			std::vector<std::shared_ptr<std::vector<BackBoneSite>>> tmploop;
			int trylen=looplen[i];
			while(tmploop.empty()){
				tmploop=BackBoneBuilder::buildlinkers(trylen,ssele[i].back(),ssele[i+1][0],std::vector<std::pair<int,int>> (),std::vector<std::pair<int,int>>(),std::set<int>());
				++maxtrytime;
				if (maxtrytime>100){
					maxtrytime=0;
					trylen++;
				}
				std::cout<<"try loop: "<<i<<" len: "<<trylen<<" time: "<<maxtrytime<<endl;
			}
			loops.push_back(tmploop);
		}
	}
	std::vector<IntrctBlck> blcks;
	int ss_count=0;
	int ss_start,ss_end;
	char ss_tmp;
	for(int i=0;i<ssele.size();i++){
		ss_start=ss_count;
		for(auto &bs:ssele[i]){
			ss_count++;
			std::string restype="GLY";
			blcks.push_back(this->mkfloatingblck(std::make_pair(-10,-10),restype,bs));
			finalss.push_back(sseq[i*2]);
			ss_tmp=sseq[i*2];
		}
		ss_end=ss_count-1;
		ss_info.push_back(std::to_string(ss_start)+"-"+std::to_string(ss_end)+" "+ss_tmp+",");
		if (ss_tmp == 'E') { ssE.push_back(std::to_string(ss_start) + "-" + std::to_string(ss_end) + ",");
		ssEH.push_back(std::to_string(ss_start) + "-" + std::to_string(ss_end) + ",");
		}
		if (ss_tmp == 'H') { ssH.push_back("0 " + std::to_string(ss_start) + " " + std::to_string(ss_end) + " 5000" + ",");
		ssEH.push_back(std::to_string(ss_start) + "-" + std::to_string(ss_end) + ",");
		}
		if (ss_tmp == 'C') { ssC.push_back(std::to_string(ss_start) + "-" + std::to_string(ss_end) + ","); }
		
	if (loop_option){
		
		if( i<ssele.size()-1){
			ss_start=ss_count;
			for(auto &bs:*(loops[i][0])){
			ss_count++;
			std::string restype="GLY";
			blcks.push_back(this->mkfloatingblck(std::make_pair(-10,-10),restype,bs));
			finalss.push_back('C');
			ss_tmp='C';
			}
		ss_end=ss_count-1;
		ss_info.push_back(std::to_string(ss_start)+"-"+std::to_string(ss_end)+" "+ss_tmp+",");
		//if (ss_tmp == 'E') { ssE.push_back(std::to_string(ss_start) + "-" + std::to_string(ss_end) + ","); }
		//if (ss_tmp == 'H') { ssH.push_back("0 "+std::to_string(ss_start) + " " + std::to_string(ss_end)+" 5000" + ","); }
		if (ss_tmp == 'C') { ssC.push_back(std::to_string(ss_start) + "-" + std::to_string(ss_end) + ","); }
		}

	}

	}
	target_->addchain(blcks);
}

IntrctMol MolModeler::splitchain(int chainid, int newcterm, int newnterm) const{
	IntrctMol imol;
	std::vector<IntrctBlck> nc1,nc2;
		for(int i=0;i<=newcterm;i++)
			nc1.push_back((*target_)({chainid,i}));
		for(int i=newnterm;i<(*target_).nresidues(chainid);i++)
			nc2.push_back((*target_)({chainid,i}));
		for(int c=0;c<(*target_).nchains();c++){
			if(c==chainid){
				if(!nc1.empty()) imol.addchain(nc1);
				if(!nc2.empty()) imol.addchain(nc2);
			} else {
				imol.addchain((*target_).getmolsystm().vec(c));
			}
	}
		imol.natoms();
		return imol;
}
IntrctBlck MolModeler::buildsidechain(const SiteIdx &site,
		const std::string &sctype,bool checkclash) {
	static PhipsiLib pplib;
	Phipsi pp = getphipsi(site);
	RotamerLib & lib = RotamerLib::rotamerlibpp(pp.phi, pp.psi);
	std::string residue = sctype;
	RotamerGroup *aagrp = lib.getAAGroup(residue);
	double emin = 10000000.0;
	Rotamer *res { nullptr };
	IntrctBlck iblksave;
	if (site.second > 0) {
		if (residue == "PRO")
			(*target_)( { site.first, site.second - 1 }).prepro = true;
		else
			(*target_)( { site.first, site.second - 1 }).prepro = false;
	}
	for (auto r : aagrp->rotList) {
		std::string rname = r->rotName;
		double rene = lib.getRotamerEnergy(rname, pplib.phipsiToIndex(&pp));
		IntrctBlck iblk = mkfloatingblck(site, r);
		if (checkclash) {
			int nc = 0;
			auto & nbors = getnearbysites(site);
			for (auto & s : nbors) {
				nc += nclashes(iblk, (*target_)(s))[1];
			}
			rene += rene + 1000 * (double) nc;
		}
		if (rene < emin) {
			emin = rene;
			iblksave = iblk;
		}
	}
	if (site.second > 0) {
		if ((*target_)(site).resname() == "PRO")
			(*target_)( { site.first, site.second - 1 }).prepro = true;
		else
			(*target_)( { site.first, site.second - 1 }).prepro = false;
	}
	assert(iblksave.topo);
	return iblksave;
}


/*
 * pocketmodeler.cpp
 *
 *  Created on: 2020年7月17日
 *      Author: hyiu
 */
#include "proteinrep/pocketmodeler.h"
#include "proteinrep/rotamerutil.h"
#include "dstl/randomengine.h"
using namespace NSPproteinrep;
#define CLASH_DISTANCE 3.5

void PocketModeler::interactivesetup(std::ostream &os, std::istream &is){
	pocket_=std::shared_ptr<IsolatedPocket>(new IsolatedPocket);
	os<< "Read some components from an existing PDB file?[Y/N]";
	std::string ans;
	getline(is, ans);
	if (ans[0] == 'Y' || ans[0] == 'y') {
		os << "Specify the input PDB filename: ";
		std::string pdbfilename;
		getline(is, pdbfilename);
		PdbReader pdb;
		try {
			pdb.readpdb(pdbfilename);
		} catch (std::exception &e) {
			std::cout << "Error occurred reading pdbfile " << pdbfilename
					<< std::endl;
			exit(0);
		}
		os <<"If all residues from the pdbfile are included, input \"ALL\""<<std::endl
				<<" Otherwise specify individual residues to include in the pocket below."
				<< std::endl
				<< "For each residue, specify its single character chainid, integer residueid"
				<< std::endl
				<< " and three-letter capital residuetype code separated by spaces in one line."
				<< std::endl << "End specification with an empty line." << std::endl;
		std::string itext;
		std::vector<IsolatedPocket::ComponentID> cmpntids;
		while (true) {
			try {
				char str[50];
				is.getline(str, 50);
				itext = std::string(str);
				if (itext.empty())
					break;
				if(itext=="ALL"){
					cmpntids.clear();
					break;
				}
				cmpntids.push_back(IsolatedPocket::ComponentID(itext));
			} catch (std::exception &e) {
				;
			}
		}
		pocket_->addcomponents(pdb, cmpntids);
	}
	std::vector<ComponentID> cmpntsfrompdb;
	for(auto &c:pocket_->components()){
		cmpntsfrompdb.push_back(c.first);
	}
	while (true) {
		os << "Add another extra residue as components?[Y/N]";
		getline(is, ans);
		if (ans[0] != 'Y' && ans[0] != 'y')
			break;
		char newchainid;
		os << "Specify its chain ID: ";
		is >> newchainid;
		getline(is, ans);
		os<< "Specify its residue type(capital three-letter) : ";
		std::string resname;
		getline(is, resname);
		AAConformer newresidue = make_aaconformer(randomrotamer(resname));
		if (pocket_->empty()) {
					pocket_->addcomponent(newresidue, newchainid,
							std::vector<IsolatedPocket::AtomID>(),
							std::vector<std::string>(), std::vector<double>());
					continue;
		}
		std::vector<IsolatedPocket::AtomID> anchors1(3);
		for (int i = 0; i < 3; i++) {
			os << "Pocket anchor atom " << i + 1 << " : ";
			is >> anchors1[i].cmpntID.chainid
					>> anchors1[i].cmpntID.residueid
					>> anchors1[i].cmpntID.resname >> anchors1[i].atmnm;
			assert(pocket_->components().find(anchors1[i].cmpntID) != pocket_->components().end());
			auto &crd=pocket_->components().at(anchors1[i].cmpntID).getglobalcrd();
			assert(crd.find(anchors1[i].atmnm) != crd.end());
			getline(is, ans);
		}
		std::vector<std::string> anchors2(3);
		std::vector<std::pair<double,double>> intcrdrngs(6);
		os << "Residue anchor atom 1: ";
		getline(is, anchors2[0]);
		os << "Distance range to pocket anchor atom 3: ";
		is >> intcrdrngs[2].first>>intcrdrngs[2].second;
		getline(is, ans);
		os << "Angle range with pocket anchor atoms 3 and 2: ";
		is >> intcrdrngs[1].first>>intcrdrngs[1].second;
		getline(is, ans);
		os << "Torsion range with  pocket anchor atoms 3, 2 and 1: ";
		is >> intcrdrngs[0].first>>intcrdrngs[0].second;
		getline(is, ans);
		os << "Residue anchor atom 2: ";
		getline(is, anchors2[1]);
		os << "Angle range with atom 1 and pocket anchor atom 3: ";
		is >> intcrdrngs[4].first>>intcrdrngs[4].second;
		getline(is, ans);
		os << "Torsion range with atom 1, pocket anchor atoms 3 and 2 ";
		is >> intcrdrngs[3].first>>intcrdrngs[3].second;
		getline(is, ans);
		os << "Residue anchor atom 3: ";
		getline(is, anchors2[2]);
		os << "Torsion range with atoms 2, 1 and pocket anchor atoms 3 and 2 ";
		is >> intcrdrngs[5].first>>intcrdrngs[5].second;
		getline(is, ans);
		double toradius = 3.14159265 / 180.0;
		for (auto & v : intcrdrngs){
			v.first *= toradius;
		   v.second *=toradius;
		}
		intcrdrngs[2].first = intcrdrngs[2].first / toradius;
		intcrdrngs[2].second=intcrdrngs[2].second/toradius;
		std::vector<double> intcrds;
		auto &rng=NSPdstl::RandomEngine<>::getinstance().realrng(0,1);
		for(auto &ic:intcrdrngs) intcrds.push_back(ic.first+rng()*(ic.second-ic.first));
		ComponentID ncid=pocket_->addcomponent(newresidue, newchainid, anchors1, anchors2,
				intcrds);
		os<<"New residues added with the following identifiers: "
				<<"chainID=\'"<<ncid.chainid<<"\'"<<" residueid= " <<ncid.residueid
				<<" residuetype= \"" << resname<<"\""<<std::endl;
		FlexRotamer frt;
		frt.cmpntid=ncid;
		frt.fxatoms=anchors2;
		flexrotamers_.push_back(frt);

		bool flexjoint=false;
		for(auto &ic:intcrdrngs){
			double diff=ic.second-ic.first;
			if(diff>1.e-3) {
				flexjoint=true;
				break;
			}
		}
		if(!flexjoint) continue;
		joints_.push_back(Joint());
		Joint &j=joints_.back();
		j.fx3atms=anchors1;
		j.mv3atms.clear();
		for(auto &a:anchors2) {
			AtomID id;
			id.cmpntID=ncid;
			id.atmnm=a;
			j.mv3atms.push_back(id);
		}
		j.mvcmpnts.insert(ncid);
		j.intcrdrngs=intcrdrngs;
	}
	chainupjoints();
	setclashexclusions();
	for(auto &c:pocket_->components()){
		cexclusions_[c.first]=std::set<ComponentID>();
	}
	for(int i=0;i<cmpntsfrompdb.size()-1;++i){
		for(int j=i+1;j<cmpntsfrompdb.size();++j){
			cexclusions_.at(cmpntsfrompdb[i]).insert(cmpntsfrompdb[j]);
			cexclusions_.at(cmpntsfrompdb[j]).insert(cmpntsfrompdb[i]);
		}
	}
}

void PocketModeler::mvaroundjoint(const Joint &joint){
	auto & rng=NSPdstl::RandomEngine<>::getinstance().realrng();
	std::vector<double> icrds;
	for(auto &ic:joint.intcrdrngs){
		icrds.push_back(ic.first+rng()*(ic.second-ic.first));
	}
	std::vector<NSPgeometry::XYZ> atms;
	for(auto &aid:joint.fx3atms) atms.push_back(
			pocket_->components()[aid.cmpntID].getglobalcrd().at(aid.atmnm));
	for(auto &aid:joint.mv3atms) atms.push_back(
				pocket_->components()[aid.cmpntID].getglobalcrd().at(aid.atmnm));
	NSPgeometry::RigidTransform rt(atms,icrds);
	for(auto &c:joint.mvcmpnts){
		auto &crds=pocket_->components()[c].getglobalcrd();
		for (auto &x:crds) rt.apply(&(x.second));
	}
}
void PocketModeler::mvrotamer(const FlexRotamer &frt){
	AAConformer *aaptr=&(pocket_->components()[frt.cmpntid]);
	AAConformer aaold=*aaptr;
	*aaptr=make_aaconformer(randomrotamer(aaold.residuename));
	aaptr->chainid_or=aaold.chainid_or;
	aaptr->residueid_or=aaold.residueid_or;
	aaptr->moveglobalcrd(aaold.getglobalcrd(),frt.fxatoms);
}

int PocketModeler::nclashes(const ComponentID &c1, const ComponentID &c2){
	auto &crd1=pocket_->components().at(c1).getglobalcrd();
	auto &crd2=pocket_->components().at(c2).getglobalcrd();
	AtomID id1,id2;
	id1.cmpntID=c1;
	id2.cmpntID=c2;
	int nc=0;
	for(auto &a1:crd1){
		id1.atmnm=a1.first;
		for(auto &a2:crd2){
			id2.atmnm=a2.first;
			double d2=(crd1.at(a1.first)-crd2.at(a2.first)).squarednorm();
			if(d2<CLASH_DISTANCE*CLASH_DISTANCE) {
				   if(clashexclusions_.at(id1).find(id2) != clashexclusions_.at(id1).end()) continue;
					nc++;
			}
		}
	}
	return nc;
}

bool PocketModeler::chainup2joints(Joint &j1, Joint &j2){
	bool mergemv{false};
	for(auto & a:j2.fx3atms){
		if(j1.mvcmpnts.find(a.cmpntID) != j1.mvcmpnts.end()) {
			mergemv=true;
			break;
		}
	}
	bool changed{false};
	if(mergemv) {
		for(auto &cm:j2.mvcmpnts){
			if(j1.mvcmpnts.find(cm) == j1.mvcmpnts.end()){
				changed=true;
				j1.mvcmpnts.insert(cm);
			}
		}
	}
	//check ring
	if(changed){
		for(auto &a:j1.fx3atms){
			if(j1.mvcmpnts.find(a.cmpntID) != j1.mvcmpnts.end()){
				std::cout<<"Ill-defined joints ( ring formed by joint-connected components) "<<std::endl;
				exit(1);
			}
		}
	}
	return changed;
}
void PocketModeler::chainupjoints(){
	bool changed=true;
	while(changed){
		changed=false;
		for(int i=0;i<joints_.size();++i){
			for(int j=0;j<joints_.size();++j){
				if(i==j) continue;
				changed=changed || chainup2joints(joints_[i],joints_[j]);
			}
		}
	}
}
static std::vector<std::string> findneighbors(std::string a1,
		const std::map<std::string, NSPgeometry::XYZ> & crds){
		std::vector<std::string> res;
		NSPgeometry::XYZ x1=crds.at(a1);
		for(auto &a2:crds){
			double r2=(x1-a2.second).squarednorm();
			if(r2<4.0)  res.push_back(a2.first);
		}
		return res;
}
void PocketModeler::setclashexclusions(){
	for(auto &c:pocket_->components()){
		for(auto &a:c.second.getglobalcrd()){
			AtomID aid;
			aid.cmpntID=c.first;
			aid.atmnm=a.first;
			clashexclusions_[aid]=std::set<AtomID>();
		}
	}
	for(auto &j:joints_){
		if(j.intcrdrngs[2].first>CLASH_DISTANCE) continue;
		std::vector<std::string> aset1=findneighbors(j.fx3atms[2].atmnm,
				pocket_->components().at(j.fx3atms[2].cmpntID).getglobalcrd());
		std::vector<std::string> aset2=findneighbors(j.mv3atms[0].atmnm,
					pocket_->components().at(j.mv3atms[0].cmpntID).getglobalcrd());
		for(auto & a1:aset1){
			AtomID id1=j.fx3atms[2];
			id1.atmnm=a1;
			for(auto &a2:aset2){
				AtomID id2=j.mv3atms[0];
				id2.atmnm=a2;
				clashexclusions_.at(id1).insert(id2);
				clashexclusions_.at(id2).insert(id1);
			}
		}
	}
}
int PocketModeler::totalclashes(){
	auto &cmp=pocket_->components();
	std::vector<ComponentID> cmpnts;
	for(auto &c:cmp) cmpnts.push_back(c.first);
	int ncl=0;
	for(int i=0;i<cmpnts.size()-1;++i){
		for(int j=i+1;j<cmpnts.size();++j){
			if(cexclusions_[cmpnts[i]].find(cmpnts[j]) != cexclusions_[cmpnts[i]].end()) continue;
			ncl+=nclashes(cmpnts[i],cmpnts[j]);
		}
	}
	return ncl;
}
bool PocketModeler::genconformation_trial(){
	for(auto &j:joints_){
		mvaroundjoint(j);
	}
	for(auto &fr:flexrotamers_){
		mvrotamer(fr);
	}
/*	ComponentID his("B 1 HIS");
	ComponentID asp("C 1 ASP");
	double d2=(pocket_->components()[his].getglobalcrd()["NE2"]-
			pocket_->components()[asp].getglobalcrd()["OD1"]).squarednorm();
	std::cout <<"Starting distance His-NE2-ASP-OD1 " << sqrt(d2)<<std::endl;*/
	int ncl=totalclashes();
	int nmv=0;
#define   MAXMOVE 1000
	while(!(ncl==0) && nmv++<MAXMOVE){
	       if(!randominternalmv(ncl)) return false;
	}
     return ncl==0;
}
bool PocketModeler::randominternalmv(int &ncl){
	int rsize=joints_.size()+flexrotamers_.size();
	auto irng=NSPdstl::RandomEngine<>::getinstance().intrng(0,rsize-1);
		int r=irng();
		if(r <joints_.size()) {
			mvaroundjoint(joints_.at(r),ncl);
		} else{
			mvrotamer(flexrotamers_.at(r-joints_.size()),ncl);
		}
	return true;
}

NSPgeometry::XYZ PocketModeler::ligandcenter() const{
	NSPgeometry::XYZ center(0,0,0);
	int nlat=0;
	if(pocket_->ligandcomponents().empty()){
		for(auto &cid: pocket_->components()){
			for(auto &crd:cid.second.getglobalcrd()){
				center =center +crd.second;
				nlat++;
			}
		}
	}else{
			for(auto &cid: pocket_->ligandcomponents()){
				for(auto &crd:pocket_->components().at(cid).getglobalcrd()){
					center =center +crd.second;
					nlat++;
				}
			}
	}
	return (1.0/(double) nlat)*center;
}
void PocketModeler::rigidtranslation(const NSPgeometry::XYZ & shift){
	for(auto &cid:pocket_->components()){
		for(auto &crd:cid.second.getglobalcrd()){
			crd.second=crd.second+shift;
		}
	}
}
void PocketModeler::rigidmv(double maxtrans,double maxrotate, const NSPgeometry::Sphere &ligandcentersphere){
	NSPgeometry::RigidTransform rt=NSPgeometry::randomrigidtransform(maxtrans,maxrotate);
	NSPgeometry::XYZ center=ligandcenter();
	assert(ligandcentersphere.insphere(center));
	while( !ligandcentersphere.insphere(rt.applytoCopy(NSPgeometry::XYZ(0,0,0))+center)){
		rt=NSPgeometry::randomrigidtransform(maxtrans,maxrotate);
	}
	this->rigidtranslation(NSPgeometry::XYZ(0.,0.,0.)-center);
	for(auto &cid:pocket_->components()){
		for(auto &crd:cid.second.getglobalcrd()){
			rt.apply(&(crd.second));
			crd.second =crd.second+center;
		}
	}
}
std::vector<NSPgeometry::XYZ> PocketModeler::copycrdto() const{
	std::vector<NSPgeometry::XYZ> res;
	for(auto &cid:pocket_->components()){
		for(auto &crd:cid.second.getglobalcrd()){
			res.push_back(crd.second);
		}
	}
	return res;
}
void PocketModeler::copycrdfrom(const std::vector<NSPgeometry::XYZ> &crds){
	int idx=0;
	for(auto &cid:pocket_->components()){
			for(auto &crd:cid.second.getglobalcrd()){
				crd.second=crds[idx++];
			}
		}
}

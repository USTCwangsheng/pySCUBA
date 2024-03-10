/*
 * intrctmol.cpp
 *
 *  Created on: 2019年12月3日
 *      Author: hyiu
 */
#include "iblock/intrctmol.h"
#include "proteinrep/aaconformer.h"
using namespace NSPintrct;
IntrctMol::IntrctMol(const IntrctMol &mola){
	*this=mola;
	for(auto &c:mol_.D_){
		for( auto &r:c){
			r.setimol(this,r.getidx2d().first,r.getidx2d().second);
		}
	}
}
void IntrctBlck::setimol(IntrctMol * imol,int id1,int id2){
	imol_=imol;
	molsystm_=&(imol->getmolsystm());
	idx2d_=std::make_pair(id1,id2);
}
void IntrctMol::setup(const std::vector<std::vector<std::string>> & sequences){
	auto &blcks=mol_.D_;
	blcks.resize(sequences.size());
	for(int c=0;c<sequences.size();++c){
		auto & chain=blcks[c];
		chain=std::vector<IntrctBlck>(sequences[c].size());
//	    std::cout <<"chain "<<c<<" size" << chain.size()<<std::endl;
		for(int r=0;r<sequences[c].size();++r){
			chain[r].setimol(this,c,r);
			this->replaceblck(NSPdstl::Idx2D{c,r},sequences[c][r]);
//			std::cout <<c<<" "<<r<<" "<<chain[r].resname()<<" "<<chain[r].natoms()<<std::endl;
		}
	}
	this->setblckoffsets(NSPdstl::Idx2D(0,0));
}
void IntrctMol::setblckoffsets(const NSPdstl::Idx2D & startblck){
	if(offsetsok_) {
#if DEBUG_LEVEL > DLEVEL_CHECKINDEX
			assert(checkblckindices());
			assert(checkaoffsets());
#endif
		return;
	}
	int aoffset=mol_[startblck].aoffset;
	residueofatoms_.resize(aoffset);
	for(int c=startblck.first;c<nchains();c++){
		int rbegin=0;
		if(c==startblck.first)rbegin=startblck.second;
		for(int r=rbegin;r<nresidues(c);++r){
			mol_(c,r).aoffset=aoffset;
			aoffset +=mol_(c,r).natoms();
			for(int i=0;i<mol_(c,r).natoms();++i)
				residueofatoms_.push_back(NSPdstl::Idx2D{c,r});
		}
	}
	offsetsok_=true;
	bsinchains=makebsinchains(*this);
	scinchains=makescinchains(*this);
	isatomfixed_.clear();
}
void IntrctMol::setup(const std::vector<std::vector<NSPproteinrep::AAConformer>> &chains,
		bool ignore_nonprotein){
	std::vector<std::vector<std::string>> sequences;
	for(auto &c:chains){
		if(ignore_nonprotein)
			if(!NSPdesignseq::ResName::resname().isStandardAminoAcid(c.at(0).residuename)) continue;
		sequences.push_back(std::vector<std::string>());
		for(auto &s:c){
			std::vector<NSPgeometry::XYZ> rcrds;
			std::string resname=s.residuename;
//			std::cout << s.residueid_or <<s.insertionid_or<<" "<<s.residuename <<std::endl;
			if(!s.sidechaincrdcomplete()) resname="GLY"; //use gly for incomplete side chains
			if(ignore_nonprotein)
				    	   if(!NSPdesignseq::ResName::resname().isStandardAminoAcid(resname)) break;
			sequences.back().push_back(resname);
		}
	}
	this->setup(sequences);
	for(int c=0;c<nchains();c++){
		for(int r=0;r<nresidues(c);++r){
			auto &blk=mol_(c,r);
			std::vector<NSPgeometry::XYZ> rcrds(blk.natoms());
			for(int a=0;a<blk.natoms();++a){
					rcrds[a]=0.1*chains[c][r].globalcrd.at(blk.atomname(a));
			}
           mol_(c,r).crds=rcrds;
		}
	}
    crds_all_.clear();
    crds_all_d_.clear();
}
/*
void IntrctMol::setup(const std::vector<std::vector<NSPproteinrep::FullSite>> &chains){
	std::vector<std::vector<std::string>> sequences;
	for(auto &c:chains){
		sequences.push_back(std::vector<std::string>());
		for(auto &s:c){
			std::vector<NSPgeometry::XYZ> rcrds;
			std::string resname=s.resname();
			if(!s.sidechaincrdcomplete()) resname="GLY"; //use gly for incomplete side chains
			sequences.back().push_back(resname);
		}
	}
	this->setup(sequences);
	for(int c=0;c<nchains();c++){
		for(int r=0;r<nresidues(c);++r){
			std::vector<NSPgeometry::XYZ> rcrds;
			if(mol_(c,r).resname()=="GLY"){
				rcrds.push_back(chains[c][r].getcrd("N"));
				rcrds.push_back(chains[c][r].getcrd("CA"));
				rcrds.push_back(chains[c][r].getcrd("C"));
				rcrds.push_back(chains[c][r].getcrd("O"));
			} else {
				chains[c][r].getcrds(rcrds);
			}
           for(auto &c:rcrds) c=0.1*c;   //Anstrom to NM
            assert(rcrds.size()==mol_(c,r).natoms());
            mol_(c,r).crds=rcrds;
		}
	}
    crds_all_.clear();
    crds_all_d_.clear();
}*/
void IntrctMol::replaceblck(NSPdstl::Idx2D idx2d,const std::string &residuename){
	IntrctBlck &blk=mol_[idx2d];
	int natoms_old=-1;
	if(blk.topo) natoms_old=blk.natoms();
	assert(NSPdesignseq::ResName::resname().isStandardAminoAcid(residuename));
	blk.topo=&(BlckTopo::getblcktopo_std(residuename));
	blk.crds.clear();
	blk.nblist_M.clear();
	blk.nblist_S.clear();
	if(natoms_old!=blk.natoms()) offsetsok_=false;
	for(auto &e:blk.energies) e=0.0;
	if(residuename=="PRO"){
		if(idx2d.second>0) mol_[idx2d].prevblck()->prepro=true;
	}
}
void IntrctMol::pasteblck(const IntrctBlck &blk){
	auto idx2d=blk.getidx2d();
	int natoms_old=(*this)(idx2d).natoms();
	std::string oldres=(*this)(idx2d).resname();
	if(blk.resname()=="PRO" && idx2d.second>0) {
		 (*this)(NSPdstl::Idx2D(idx2d.first,idx2d.second-1)).prepro=true;
    }  else if(oldres=="PRO" && idx2d.second>0 ){
	    (*this)(NSPdstl::Idx2D(idx2d.first,idx2d.second-1)).prepro=false;
	}
	int offset_old=(*this)(idx2d).aoffset;
	(*this)(idx2d)=blk;
	(*this)(blk.getidx2d()).aoffset=offset_old;
	if(natoms_old != blk.natoms()) {
		offsetsok_=false;
	}
	crds_all_.clear();
	crds_all_d_.clear();
}
void IntrctMol::changecrds(const std::vector<NSPgeometry::XYZ> &newcrds) const{
#if DEBUG_LEVEL > DLEVEL_CHECKINDEX
			assert(checkblckindices());
			assert(checkaoffsets());
#endif
	assert(offsetsok_);
	assert(newcrds.size()>=natoms());
	crds_all_=newcrds;
	crds_all_d_=NSPgeometry::XYZvtodoublev(crds_all_);
	for(int  c=0;c<nchains();++c){
		for(int r=0;r<nresidues(c);++r){
			const IntrctBlck &b=mol_(c,r);
			b.crds.resize(b.natoms());
			std::copy(newcrds.begin()+b.aoffset,newcrds.begin()+b.aoffset+b.natoms(),
					b.crds.begin());
		}
	}
}
const std::vector<NSPgeometry::XYZ> & IntrctMol::recollectcrds_all() const {
#if DEBUG_LEVEL > DLEVEL_CHECKINDEX
			assert(checkblckindices());
			assert(checkaoffsets());
#endif
	assert(offsetsok_);
	if(crds_all_.size()==natoms()) {
#if DEBUG_LEVEL > DLEVEL_CHECKCRDS
			assert(checkcrds());
#endif
		return crds_all_;
	}
	crds_all_.resize(this->natoms());
	for(int  c=0;c<nchains();++c){
		for(int r=0;r<nresidues(c);++r){
			const IntrctBlck &b=mol_(c,r);
			std::copy(b.crds.begin(),b.crds.end(),crds_all_.begin()+b.aoffset);
		}
	}
	crds_all_d_=NSPgeometry::XYZvtodoublev(crds_all_);
	return crds_all_;
}
void IntrctMol::translatechaincenter( const NSPgeometry::XYZ & newcenter, int chainid){
	NSPgeometry::XYZ oldcenter(0.0,0.0,0.0);
	int nat=0;
	for(int r=0; r<nresidues(chainid); ++r){
		for(int a=0; a<mol_(chainid,r).natoms();++a){
			oldcenter =oldcenter+mol_(chainid,r).crds[a];
			++nat;
		}
	}
	oldcenter=oldcenter*(1.0/(double) nat);
	NSPgeometry::XYZ t=newcenter-oldcenter;
	for(int r=0; r<nresidues(chainid); ++r){
		for(int a=0; a<mol_(chainid,r).natoms();++a){
			mol_(chainid,r).crds[a] = mol_(chainid,r).crds[a] +t;
		}
	}
	crds_all_.clear();
}
void IntrctMol::clear_energies(){
	for(int  c=0;c<nchains();++c){
		for(int r=0;r<nresidues(c);++r){
			for(auto &e:mol_(c,r).energies)e=0.0;
		}
	}
}
void IntrctMol::mksteric_nblists(const IntrctPara &param) const{
		for(auto & c:mol_.D_){
			for(auto &r:c){
				r.mksteric_nblist(param);
			}
		}
}
void IntrctMol::forces_all(const IntrctPara &param,
		std::vector<NSPgeometry::XYZ> &dedx) const{
	recollectcrds_all();
	this->calcphipsicodes(param.phipsicodes_updateonly);
	if(param.calc_sscodes) this->calcsscodes(param.sscodes_updateonly);
	if(param.calc_nblists) this->mksteric_nblists(param);
	if(param.enedetails){
		auto & oss=Results::ofstreams(param.jobname,"SecondaryStr");
		auto ssstrings=this->sscodestrings();
		for(auto &s:ssstrings) oss<<s<<std::endl;
	}
	this->calcconformercodes();
	dedx.assign(natoms(),NSPgeometry::XYZ{0,0,0});
	if(param.weight_coval >0) forces_cov(param,dedx);
    if(param.weight_steric>0) forces_steric(param,dedx);
	if(param.weight_phipsi>0) forces_phipsi(param,dedx);
	if(param.weight_localstr>0) forces_localstr(param,dedx);
	if(param.weight_localhb>0) forces_localhb(param,dedx);
	if(param.weight_sitepair>0) forces_sitepair(param,dedx);
	if(param.weight_rotamer>0) forces_rotamer(param,dedx);
	if(param.weight_scpacking>0) forces_scpacking(param,dedx);
}
double IntrctMol::sum_energies(std::array<double,IntrctBlck::ENESIZE> &energies) const {
	for(auto &e:energies) e=0.0;
	for(int  c=0;c<nchains();++c){
		for(int r=0;r<nresidues(c);++r){
			for(int i=0;i<IntrctBlck::ENESIZE;++i){
				energies[i]+=mol_(c,r).energies[i];
			}
		}
	}
	double esum=0.0;
	for(auto &e:energies)
		esum +=e;
	return esum;
}
const std::vector<bool> & IntrctMol::setatomfixed() {
	if(!offsetsok_) {
		setblckoffsets(NSPdstl::Idx2D{0,0});
	}
	if(isatomfixed_.size()==natoms()) {
#if DEBUG_LEVEL > DLEVEL_CHECKINDEX
			assert(checkatomsfixed());
#endif
		return isatomfixed_;
	}
	isatomfixed_.assign(natoms(),false);
	auto & ratoms=residueofatoms();
	for(int i=0;i<natoms();++i){
		auto &idx2d=ratoms[i];
		auto &blk=(*this)(idx2d);
		if(blk.activemod==IntrctBlck::ALL) continue;
		if(blk.activemod==IntrctBlck::INACTIVE) {
			isatomfixed_[i]=true;
		}
		else if(blk.activemod==IntrctBlck::SIDECHAIN){
			if(blk.topo->atomips[i-blk.aoffset].ismainchain) isatomfixed_[i]=true;
		}
	}
	return isatomfixed_;
}
std::vector<int> IntrctMol::atomindices(const std::map<int,std::set<int>> & chainblcks) const{
		   std::vector<int> aindices;
		   assert(offsetsok_);
		   for(auto &cb:chainblcks){
			   for(auto bidx:cb.second){
				   auto &blck=this->getblck(NSPdstl::Idx2D{cb.first,bidx});
				   int aoffset=blck.aoffset;
				   for(int i=0;i<blck.natoms();++i) aindices.push_back(aoffset+i);
			   }
		   }
		   return aindices;
}
void IntrctMol::writepdb(std::ofstream &ofs) const{
	for(int ic=0;ic<nchains();++ic){
		for(int r=0;r<nresidues(ic);++r){
			std::vector<NSPproteinrep::PdbRecord> records=
					(*this)(NSPdstl::Idx2D{ic,r}).topdbrecords();
			for(auto &rc:records)
				ofs<<rc.toString()<<std::endl;
		}
	}
}
#include "iblock/molmodeler.h"
void IntrctMol::ssregions(int chainid, std::vector<std::pair<int,int>> &helixregions,
		 std::vector<std::pair<int,int>> &strandregions){
		std::vector<char> ssseq;
		if(this->sscodes.empty()) this->calcsscodes();
		for(auto &code:sscodes.at(chainid)){
			ssseq.push_back(code.charcode());
		}
		NSPintrct::ssregions(ssseq,helixregions,strandregions);
}

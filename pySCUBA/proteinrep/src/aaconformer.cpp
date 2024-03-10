/*
 * aaconformer.cpp
 *
 *  Created on: 2017年6月27日
 *      Author: hyliu
 */

#include "proteinrep/aaconformer.h"
#include "proteinrep/rotamerutil.h"
#include <set>
#include <cassert>
using namespace NSPproteinrep;
using namespace NSPgeometry;
std::vector<std::string> AAConformer::mainchainatoms{"N","CA","C","O"};
std::vector<std::string> AAConformer::backboneatoms{"N","CA","C"};
std::map<std::string,std::vector<std::string>> AAConformer::sidechainatoms{
	{"GLY",{}},
	{"ALA",{"CB"}},
	{"CYS",{"CB","SG"}},
	{"ASP",{"CB","CG","OD1","OD2"}},
	{"GLU",{"CB","CG","CD","OE1","OE2"}},
	{"PHE",{"CB","CG","CD1","CD2","CE1","CE2","CZ"}},
	{"HIS",{"CB","CG","ND1","CD2","CE1","NE2"}},
	{"ILE",{"CB","CG1","CG2","CD1"}},
	{"LYS",{"CB","CG","CD","CE","NZ"}},
	{"LEU",{"CB","CG","CD1","CD2"}},
	{"MET",{"CB","CG","SD","CE"}},
	{"ASN",{"CB","CG","OD1","ND2"}},
	{"PRO",{"CB","CG","CD"}},
	{"GLN",{"CB","CG","CD","OE1","NE2"}},
	{"ARG",{"CB","CG","CD","NE","CZ","NH1","NH2"}},
	{"SER",{"CB","OG"}},
	{"THR",{"CB","OG1","CG2"}},
	{"VAL",{"CB","CG1","CG2"}},
	{"TRP",{"CB","CG","CD1","CD2","NE1","CE2","CE3","CZ2","CZ3","CH2"}},
	{"TYR",{"CB","CG","CD1","CD2","CE1","CE2","CZ","OH"}}
};
std::map<std::string,NSPgeometry::XYZ> AAConformer::moveglobalcrd(
		        const std::map<std::string,NSPgeometry::XYZ> &refgcrd,
	    const std::vector<std::string> &frmatoms){
        calclocalcrd(frmatoms);
		localframe=calclocalframe(refgcrd,frmatoms);
		return calcglobalcrd();
}
NSPgeometry::LocalFrame AAConformer::calclocalframe(
		const std::map<std::string, NSPgeometry::XYZ> & gcrd,
		std::vector<std::string> frmatoms) {
	if(frmatoms.empty()){
		frmatoms.push_back("CA");
		frmatoms.push_back("N");
		frmatoms.push_back("C");
	}
	assert(frmatoms.size()==3);
	XYZ origin = gcrd.at(frmatoms[0]);
	XYZ pn = gcrd.at(frmatoms[1]);
	XYZ pc = gcrd.at(frmatoms[2]);
	XYZ pcb;
//	if(gcrd.find("CB") != gcrd.end())pcb=gcrd.at("CB");
	pcb=NSPgeometry::InternaltoXYZ(origin,pc,pn,1.5,109.5*3.14159265/180.0,
			120*3.14159265/180.0);
	XYZ pxy=0.5*(pn+pc);
	return make_localframe(origin, pcb, pxy);
}

std::map<std::string, NSPgeometry::XYZ> & AAConformer::calclocalcrd(
		const std::vector<std::string> &frameatoms) {
	localframe = calclocalframe(globalcrd,frameatoms);
	localcrd.clear();
	for (auto &p : globalcrd) {
		XYZ lx = localframe.global2localcrd(p.second);
		localcrd.insert(std::make_pair(p.first, lx));
	}
	return localcrd;
}
std::map<std::string, NSPgeometry::XYZ> & AAConformer::calcglobalcrd() {
	globalcrd.clear();
	for (auto &p : localcrd) {
		XYZ gx = localframe.local2globalcrd(p.second);
		globalcrd.insert(std::make_pair(p.first, gx));
	}
	return globalcrd;
}
bool NSPproteinrep::ismainchain(const std::string & atomname){
	for(auto & m:AAConformer::mainchainatoms) if(atomname == m) return true;
	if(atomname=="H") return true;
	return false;
}
bool NSPproteinrep::isbackbone(const std::string & atomname){
	for(auto & b:AAConformer::backboneatoms) if(atomname == b) return true;
	return false;
}
double NSPproteinrep::sidechainRMSD2_local(const AAConformer &c1,
		const AAConformer &c2) {
	assert(c1.residuename == c2.residuename);
	const std::map<std::string, XYZ> & lcrd1 = c1.getlocalcrd();
	const std::map<std::string, XYZ> & lcrd2 = c2.getlocalcrd();
	double res = 0.0;
	int nside = 0;
	for (auto &p1 : lcrd1) {
		if (ismainchain(p1.first))
			continue;
		++nside;
		res += distance2(p1.second, lcrd2.at(p1.first));
	}
	return res / (double) nside;
}

double NSPproteinrep::RMSD2_global(const AAConformer &c1, const AAConformer &c2) {
	assert(c1.residuename == c2.residuename);
	const std::map<std::string, XYZ> & gcrd1 = c1.getglobalcrd();
	const std::map<std::string, XYZ> & gcrd2 = c2.getglobalcrd();
	double res = 0.0;
	for (auto &p1 : gcrd1) {
		res += distance2(p1.second, gcrd2.at(p1.first));
	}
	return res / (double) gcrd1.size();
}
double  NSPproteinrep::backboneRMSD2(const AAConformer &c1, const AAConformer &c2) {
	assert(c1.residuename == c2.residuename);
	const std::map<std::string, XYZ> & gcrd1 = c1.getglobalcrd();
	const std::map<std::string, XYZ> & gcrd2 = c2.getglobalcrd();
	double res = 0.0;
	int nback = 0;
	for (auto &p1 : gcrd1) {
		if (!isbackbone(p1.first))
			continue;
		++nback;
		res += distance2(p1.second, gcrd2.at(p1.first));
	}
	return res / (double) nback;
}
AAConformer AAConformer::rotamerconformer() const{
	BackBoneSite bs=make_backbonesite();
	auto &rotamers=RotamerDistribution::getdistribution(residuename).rotamers;
	AAConformer res;
	double rmsd2min=10000.0;
	for(auto r:rotamers){
		AAConformer v2=make_aaconformer(r,&bs);
		double rmsd2=RMSD2_global(*this,v2);
		if(rmsd2<rmsd2min){
			res=v2;
			rmsd2min=rmsd2;
		}
	}
    return res;
}
AAConformer NSPproteinrep::make_aaconformer(const std::vector<PdbRecord> & residuerecords,
		std::vector<AAConformer> *altconformers, bool replaceMSE){
	AAConformer aa;
	aa.residuename=residuerecords[0].residuename;
	bool mse2met=replaceMSE && aa.residuename=="MSE";
	if(mse2met) aa.residuename="MET";
	aa.chainid_or=residuerecords[0].chainid;
	aa.residueid_or=residuerecords[0].residueid;
	aa.insertionid_or=residuerecords[0].insertionid;
	std::set<char> states;
	for(auto & a:residuerecords){
		if(a.conformerid != ' ')
		states.insert(a.conformerid);
	}
	if(states.empty()) states.insert(' ');
	std::vector<char> states_v;
	for(auto s:states) states_v.push_back(s);
	if(altconformers && states_v.size()>1) {
		altconformers->clear();
		for(int i=0; i<states_v.size()-1;++i) {
			altconformers->push_back(AAConformer());
			altconformers->back().residuename=aa.residuename;
		}
	}
	for(auto & a:residuerecords) {
		if(a.conformerid == states_v[0] ||a.conformerid==' ') {
//			assert(aa.globalcrd.find(a.atomname) == aa.globalcrd.end());
			if(aa.globalcrd.find(a.atomname) == aa.globalcrd.end()){
				aa.globalcrd.insert(std::make_pair(a.atomname,XYZ(a.x,a.y,a.z)));
				if(mse2met && a.atomname=="SE")  aa.atomlist.push_back("SD");
				else aa.atomlist.push_back(a.atomname);
			}
		}
		else if(altconformers) {
			int id;
			for(int i=1; i<states_v.size();++i) if(a.conformerid == states_v[i]) id=i-1;
			altconformers->at(id).globalcrd.insert(std::make_pair(a.atomname,XYZ(a.x,a.y,a.z)));
		}
	}
	return aa;
}
std::vector<PdbRecord> AAConformer::make_pdbrecords(char chainid,int resid,char insertionid,int atomid0) const {
	std::vector<PdbRecord> res;
	std::vector<std::string> atomnames;
	if(isaa()){
		for(auto &nm:mainchainatoms) atomnames.push_back(nm);
		if(isnaturalaa()){
			for(auto &nm:sidechainatoms.at(residuename))
				atomnames.push_back(nm);
		}
	} else {
		for(auto & atm:globalcrd) atomnames.push_back(atm.first);
	}
	for(auto & atm:atomnames) {
		if(globalcrd.find(atm) == globalcrd.end()) continue;
		res.push_back(PdbRecord());
		PdbRecord & rec=res.back();
		if(isnaturalaa())rec.label="ATOM";
		else rec.label="HETATM";
		rec.chainid=chainid;
		rec.atomname=atm;
		rec.atomid=atomid0++;
		rec.residuename=residuename;
		rec.namesymbol=atm.substr(0,1);
		rec.elementname[0]=' ';
		rec.elementname[1]=atm[0];
		rec.namemodifier=atm.substr(1);
		rec.residueid=resid;
		rec.insertionid=insertionid;
		rec.x=globalcrd.at(atm).x_;
		rec.y=globalcrd.at(atm).y_;
		rec.z=globalcrd.at(atm).z_;
	}
	return res;
}
BackBoneSite AAConformer::make_backbonesite() const {
	BackBoneSite bs;
	bs.resid=this->residueid_or;
	bs.resseq=this->residueid_or;
	bs.resname=this->residuename;
	bs.data_[BackBoneSite::PHI] = 360.0;
	double p=NSPgeometry::torsion(globalcrd.at("N"), globalcrd.at("CA"),
			globalcrd.at("C"), globalcrd.at("O"))*180.0/3.14159265+180.0;
	if(p>180.0) p -=360.0;
	bs.data_[BackBoneSite::PSI] =p;
	bs.data_[BackBoneSite::OMIGA] = 180.0;
	bs.changecrd(this->getbackbonecrd());
	return bs;
}
int AAConformer::removeatomsnotlisted(){
	assert(isnaturalaa());
	std::vector<std::string> removes;
	for(auto &atm:globalcrd){
		std::string name=atm.first;
		bool toremove=true;
		for(auto &n:mainchainatoms) if(n==name){toremove=false;break;}
		if(toremove)
			for(auto & n:sidechainatoms[residuename]) if(n==name) {toremove=false;break;}
		if(toremove) removes.push_back(name);
	}
	for(auto & n:removes) {globalcrd.erase(n); localcrd.erase(n);}
	return removes.size();
}
AAConformer NSPproteinrep::make_aaconformer(NSPdesignseq::Rotamer *rt,
		const BackBoneSite *site){
	   BackBoneSite bs;
	   if(!site) {
		   genbackbonesite(nullptr, false,180.0, 120.0, &bs);
	   } else {
		   bs=*site;
	   }
	   auto cs=NSPdesignseq::getBackboneSiteLocalFrame(bs);
	   std::vector<NSPgeometry::XYZ> sccrds;
	   	rt->buildSidechain(cs, sccrds);
	   	AAConformer conf;
	   	conf.residuename=rt->triName;
	   	conf.chainid_or='A';
	   	conf.residueid_or=1;
	   	conf.globalcrd["N"]=bs.ncrd();
	   	conf.atomlist.push_back("N");
	   	conf.globalcrd["CA"]=bs.cacrd();
	 	conf.atomlist.push_back("CA");
	   	conf.globalcrd["C"]=bs.ccrd();
	 	conf.atomlist.push_back("C");
	   	conf.globalcrd["O"]=bs.ocrd();
	 	conf.atomlist.push_back("O");
	   	for(int i=0;i<rt->atomNameList.size();++i){
	   		conf.globalcrd[rt->atomNameList[i]]=sccrds[i];
	   		conf.atomlist.push_back(rt->atomNameList[i]);
	   	}
	   return conf;
}

void AAConformersInModel::readpdbfile(const std::string &filename){
	PdbReader pdb;
	pdb.readpdb(filename);
	getconformers(pdb);
}
void AAConformersInModel::getconformers(const PdbReader &pdb){
	mappdbkeyint=pdb.mappdbkeyint();
	conformers.clear();
	altconformers.clear();
	int chainnumber=0;
	for(auto & c:pdb.records()) {
		assert(chainnumber++ == mappdbkeyint->chainNumber(c.first));
		conformers.push_back(std::vector<AAConformer>());
		altconformers.push_back(std::vector<std::vector<AAConformer>>());
		std::vector<AAConformer> & residues=conformers.back();
		std::vector<std::vector<AAConformer>> & alt=altconformers.back();
		int posi=0;
		for(auto &r:c.second) {
//			if(r.second[0].label != "ATOM") continue;
//			std::cout <<"---"<<r.second[0].residueid <<r.second[0].insertionid<<" "<<r.second[0].residuename<<std::endl;
			alt.push_back(std::vector<AAConformer>());
			residues.push_back(make_aaconformer(r.second,&(alt.back())));
			assert(posi++ ==
					mappdbkeyint->posiNumber({residues.back().residueid_or,residues.back().insertionid_or},c.first));
		}
	}
}

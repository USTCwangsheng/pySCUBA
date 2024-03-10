/*
 * blcktopo.cpp
 *
 *  Created on: 2019年12月1日
 *      Author: hyiu
 */

#include "iblock/vsctypes.h"
#include "iblock/backboneff.h"
#include "iblock/blcktopo.h"

using namespace NSPintrct;
#include "iblock/vsctypes.h"
#include "iblock/blcktopo.h"
#include "iblock/internalcrds.h"
#include "geometry/calculators.h"

using namespace NSPintrct;
using namespace NSPproteinrep;
using namespace NSPgeometry;

BlckTopo mkblcktopo_std(const std::string &residuename) {
	BlckTopo topo;
        BackBoneFF backboneff;
	const VSCType & vsc = VSCType::getVSCType(residuename);
	topo.blckname = residuename;
	topo.natoms = vsc.nscatoms + 4;
	topo.atomnames.push_back("N");
	topo.atomnames.push_back("CA");
	for (auto &a : vsc.atomnames)
		topo.atomnames.push_back(a);
	topo.atomnames.push_back("C");
	topo.atomnames.push_back("O");
	for(int i=0;i<topo.atomnames.size();++i){
		topo.atomseqs[topo.atomnames[i]]=i;
	}
	std::vector<AtomIdx3D> idx3ds;
	for (int i = 0; i < topo.natoms; ++i) {
		idx3ds.push_back(AtomIdx3D { 0, 0, i });
	}
	std::array<AtomIdx3D, 2> bidx;
	bidx[0] = idx3ds[0];
	bidx[1] = idx3ds[1];
	topo.bonds.push_back(
			Bond(bidx, backboneff.b0_nca*A2NM, 2 * KBT * backboneff.kb_nca/(A2NM*A2NM)));
	for (int b = 0; b < vsc.newbonds.size(); ++b) {
		bidx[0] = idx3ds[vsc.newbonds[b].first];
		bidx[1] = idx3ds[vsc.newbonds[b].second];
		topo.bonds.push_back(Bond(bidx, vsc.b0[b]*A2NM, vsc.kb0[b]/(A2NM*A2NM)));
	}
	bidx[0] = idx3ds[1];
	bidx[1] = idx3ds[topo.natoms - 2];
	topo.bonds.push_back(
			Bond(bidx, backboneff.b0_cac*A2NM, 2 * KBT * backboneff.kb_cac/(A2NM*A2NM)));
	bidx[0] = idx3ds[topo.natoms - 2];
	bidx[1] = idx3ds[topo.natoms - 1];
	topo.bonds.push_back(
			Bond(bidx, backboneff.b0_co*A2NM, 2 * KBT * backboneff.kb_co/(A2NM*A2NM)));
	//bond to next residue
	bidx[0] = idx3ds[topo.natoms - 2];
	bidx[1] = AtomIdx3D { 0, 1, 0 };
	topo.bonds.push_back(
			Bond(bidx, backboneff.b0_cn*A2NM, 2 * KBT * backboneff.kb_cn/(A2NM*A2NM)));

	std::array<AtomIdx3D, 3> aidx;
	aidx[0] = idx3ds[0];
	aidx[1] = idx3ds[1];
	aidx[2] = idx3ds[topo.natoms - 2];
	topo.angles.push_back(
			Angle(aidx, cos(backboneff.t0_ncac*3.14159265/180.0),
					2 * KBT * KANG_FAC * backboneff.kt_ncac));
	aidx[0] = idx3ds[1];
	aidx[1] = idx3ds[topo.natoms - 2];
	aidx[2] = idx3ds[topo.natoms - 1];
	topo.angles.push_back(
			Angle(aidx, cos(backboneff.t0_caco*3.14159265/180.0),
					2 * KBT * KANG_FAC * backboneff.kt_caco));
	for (int a = 0; a < vsc.newangles.size(); ++a) {
		for (int i = 0; i < 3; ++i)
			aidx[i] = idx3ds[vsc.newangles[a][i]];
		topo.angles.push_back(Angle(aidx, cos(vsc.a0[a]*3.14159265/180.0), vsc.ka0[a]));
	}
	aidx[0] = idx3ds[1];
	aidx[1] = idx3ds[topo.natoms - 2];
	aidx[2] = AtomIdx3D { 0, 1, 0 };
	topo.angles.push_back(
			Angle(aidx, cos(backboneff.t0_cacn*3.14159265/180.0),
					2 * KBT * KANG_FAC * backboneff.kt_cacn));
	aidx[0] = idx3ds[topo.natoms - 1];
	aidx[1] = idx3ds[topo.natoms - 2];
	aidx[2] = AtomIdx3D { 0, 1, 0 };
	topo.angles.push_back(
			Angle(aidx, cos(backboneff.t0_ocn*3.14159265/180.0),
					2 * KBT * KANG_FAC * backboneff.kt_ocn));
	aidx[0]=idx3ds[topo.natoms - 2];
	aidx[1] = AtomIdx3D { 0, 1, 0 };
	aidx[2]=AtomIdx3D {0, 1, 1};
	topo.angles.push_back(
			Angle(aidx, cos(backboneff.t0_cnca*3.14159265/180.0),
					2 * KBT * KANG_FAC * backboneff.kt_cnca));
	std::array<AtomIdx3D, 4> pidx;
	for (int a = 0; a < vsc.newimpdihs.size(); ++a) {
		for (int i = 0; i < 4; ++i)
			pidx[i] = idx3ds[vsc.newimpdihs[a][i]];
		topo.impdihs.push_back(ImpDih(pidx, vsc.imp0[a]*3.14159265/180.0, vsc.kimp0[a]));
	}
    pidx[0] = idx3ds[topo.natoms - 2];
	pidx[1] = AtomIdx3D { 0, 1, 0 };
	pidx[2] = idx3ds[1];
	pidx[3] = idx3ds[topo.natoms - 1];
	topo.impdihs.push_back(
			ImpDih(pidx, backboneff.p_cncao*3.14159265/180.0,
					2 * KBT * KANG_FAC * backboneff.kp_cncao));
	pidx[0] = idx3ds[1];
	pidx[1] = idx3ds[topo.natoms - 2];
	pidx[2] = AtomIdx3D { 0, 1, 0 };
	pidx[3] = AtomIdx3D { 0, 1, 1 };
	topo.pepimpdih1=topo.impdihs.size();
	topo.impdihs.push_back(
			ImpDih(pidx, backboneff.p_cacnca*3.14159265/180.0,
					2 * KBT * KANG_FAC * backboneff.kp_cacnca));
	pidx[0] = idx3ds[topo.natoms - 1];
	topo.pepimpdih2=topo.impdihs.size();
	topo.impdihs.push_back(
			ImpDih(pidx, backboneff.p_ocnca*3.14159265/180.0,
					2 * KBT * KANG_FAC * backboneff.kp_ocnca));
	pidx[0] = idx3ds[0];
	pidx[1] = idx3ds[1];
	pidx[2] = idx3ds[topo.natoms - 2];
	pidx[3] = AtomIdx3D { 0, 1, 0 };
	topo.torsions.push_back(Torsion { pidx[0], pidx[1], pidx[2], pidx[3] });
	pidx[0] = idx3ds[0];
	pidx[1] = idx3ds[1];
	pidx[2] = idx3ds[topo.natoms - 2];
	pidx[3] = idx3ds[topo.natoms - 1];
	topo.torsions.push_back(Torsion { pidx[0], pidx[1], pidx[2], pidx[3] });
	pidx[0] = idx3ds[topo.natoms - 2];
	pidx[1] = AtomIdx3D { 0, 1, 0 };
	pidx[2] = AtomIdx3D { 0, 1, 1 };
	pidx[3] = AtomIdx3D { 0, 1, -2 }; //-2 represents C
	topo.torsions.push_back(Torsion { pidx[0], pidx[1], pidx[2], pidx[3] });
	pidx[1] = AtomIdx3D { 0, 1, 0 };
	pidx[2] = AtomIdx3D { 0, 1, 1 };
	pidx[3] = AtomIdx3D { 0, 1, -3 };  //-3 represents CB
	topo.torsions.push_back(Torsion { pidx[0], pidx[1], pidx[2], pidx[3] });
	for (int a = 0; a < vsc.newtorsions.size(); ++a) {
		bool keep = true;
		int idx = 0;
		for (auto k : vsc.newtorsions[a]) {
			if (k < 0) {
				keep = false;
				break;
			}
			if (k < topo.natoms) {
				pidx[idx++] = idx3ds[k];
			} else {
				pidx[idx++] = AtomIdx3D { 0, 1, k - topo.natoms };
			}
		}
		if (!keep)
			continue;
		topo.torsions.push_back(Torsion { pidx[0], pidx[1], pidx[2], pidx[3]});
	}
	topo.atomips.push_back(
			AtomIP(VSCType::getstericatomtype(vsc.resname, "N"),
					backboneff.sigma_n*A2NM, backboneff.eps, backboneff.sigma_nhb*A2NM,
					true, false, true, false, true));
	topo.atomips.push_back(
			AtomIP(VSCType::getstericatomtype(vsc.resname, "CA"),
					backboneff.sigma_ca*A2NM, backboneff.eps, backboneff.sigma_ca*A2NM,
					false, false, true, false, true));
	for (auto & a : vsc.atomnames) {
		int atype = VSCType::getstericatomtype(vsc.resname, a);
		double sigma = VSCType::packingatomtypes[atype].radius;
		int hbtype = VSCType::packingatomtypes[atype].hbtype;
		bool hbdonor = (hbtype == 1 || hbtype == 3);
		bool hbacceptor = (hbtype == 2 || hbtype == 3);
		double sigmahb = 2.8;
		topo.atomips.push_back(
				AtomIP(atype, sigma*A2NM, 0.2, sigmahb*A2NM, hbdonor, hbacceptor, false,
						true, true));
	}
	topo.atomips.push_back(
			AtomIP(VSCType::getstericatomtype(vsc.resname, "C"),
					backboneff.sigma_c*A2NM, backboneff.eps, backboneff.sigma_c*A2NM,
					false, false, true, false, true));
	topo.atomips.push_back(
			AtomIP(VSCType::getstericatomtype(vsc.resname, "O"),
					backboneff.sigma_o*A2NM, backboneff.eps, backboneff.sigma_ohb*A2NM,
					false, true, true, false, true));


	topo.excld.assign(topo.natoms,std::set<AtomIdx3D>());
	for(auto & b:topo.bonds) setexcld<2>(topo.excld,b.aidx);
	for(auto & a:topo.angles) setexcld<3>(topo.excld,a.aidx);
	for(auto & t:topo.impdihs) setexcld<4>(topo.excld,t.aidx);
	for(auto & t:topo.torsions) setexcld<4>(topo.excld,t);
	topo.nblist_internal_M.assign(topo.natoms,std::set<int>());
	topo.nblist_internal_S.assign(topo.natoms,std::set<int>());
	for(int i=0;i<topo.natoms-1;++i){
		for(int j=i+1;j<topo.natoms;++j){
				if(topo.excld[i].find(idx3ds[j])==topo.excld[i].end()){
					if(topo.atomips[i].issidechain||topo.atomips[j].issidechain)
						topo.nblist_internal_S[i].insert(j);
					else
						topo.nblist_internal_M[i].insert(j);
				}
		}
	}
	return topo;
}
const BlckTopo & BlckTopo::getblcktopo_std(const std::string &residuename){
	static std::map<std::string, BlckTopo> residuetopos;
	if(residuetopos.find(residuename)== residuetopos.end())
		residuetopos[residuename]=mkblcktopo_std(residuename);
	return residuetopos.at(residuename);
}
const BlckTopo &BlckTopo::getblcktopo_nonstd(const std::string &ligandname,
		const NSPproteinrep::AAConformer *ligandconformer){
	static std::map<std::string,BlckTopo> ligandtopos;
	if(ligandtopos.find(ligandname) != ligandtopos.end()) return ligandtopos.at(ligandname);
	assert(ligandconformer != nullptr);
	ligandtopos[ligandname]=BlckTopo();
	auto &topo=ligandtopos[ligandname];
	std::vector<std::string> atomnames;
	std::vector<XYZ> crds;
	for(auto &c:ligandconformer->getglobalcrd()) {
		atomnames.push_back(c.first);
		crds.push_back(c.second);
	}
	topo.atomnames=atomnames;
	topo.natoms=atomnames.size();
	for(int i=0;i<topo.natoms;++i){
		int atype=VSCType::getstericatomtype("LIG", "ANY");
		if(atype<0){
			std::cout <<"Cannot assign ligand atom type. Make sure you have the \"ligandatomtypes.dat\"\n"
					<< "file in your working directory or in the SCUBA_DATAPATH directory\n";
			exit(1);
		}
		double sigma = VSCType::packingatomtypes[atype].radius;
		int hbtype = VSCType::packingatomtypes[atype].hbtype;
		bool hbdonor = (hbtype == 1 || hbtype == 3);
		bool hbacceptor = (hbtype == 2 || hbtype == 3);
		double sigmahb = 2.8;
		topo.atomips.push_back(
				AtomIP(atype, sigma*A2NM, 0.2, sigmahb*A2NM, hbdonor, hbacceptor, false,
						false, false));
	}
	std::vector<AtomIdx3D> idx3ds;
	for (int i = 0; i < topo.natoms; ++i) {
		idx3ds.push_back(AtomIdx3D { 0, 0, i });
	}
	InternalCrds ic(crds);
	std::array<AtomIdx3D, 2> bidx;
	double kb=2.0*KBT*1000.0/(A2NM*A2NM);
	for(auto &b:ic.bonds()){
			bidx[0]=idx3ds[b.first];
			bidx[1]=idx3ds[b.second];
			double b0=distance(crds[b.first],crds[b.second]);
			topo.bonds.push_back(Bond(bidx, b0,kb));
	}
	std::array<AtomIdx3D,3> aidx;
	double ka=2 * KBT * KANG_FAC*0.05;
	for(auto &a:ic.angles()){
		for(int i=0;i<3;++i) aidx[i]=idx3ds[a[i]];
		double a0=angle(crds[a[0]],crds[a[1]],crds[a[2]]);
		topo.angles.push_back(Angle(aidx,cos(a0),ka));
	}
	std::array<AtomIdx3D,4> iidx;
	double ki=2 * KBT * KANG_FAC*0.01;
	for(auto &p:ic.impdihs()){
		for(int i=0;i<4;++i) iidx[i]=idx3ds[p[i]];
		double imp0=torsion(crds[p[0]],crds[p[1]],crds[p[2]],crds[p[3]]);
		topo.impdihs.push_back(ImpDih(iidx,imp0,ki));
	}
	for(auto &p:ic.torsions()){
		for(int i=0;i<4;++i) iidx[i]=idx3ds[p[i]];
		topo.torsions.push_back(Torsion{iidx[0],iidx[1],iidx[2],iidx[3]});
	}
	topo.excld.assign(topo.natoms,std::set<AtomIdx3D>());
	for(auto & b:topo.bonds) setexcld<2>(topo.excld,b.aidx);
	for(auto & a:topo.angles) setexcld<3>(topo.excld,a.aidx);
	for(auto & t:topo.impdihs) setexcld<4>(topo.excld,t.aidx);
	for(auto & t:topo.torsions) setexcld<4>(topo.excld,t);
	topo.nblist_internal_M.assign(topo.natoms,std::set<int>());
	topo.nblist_internal_S.assign(topo.natoms,std::set<int>());
	ligandtopos[ligandname]=topo;
	return ligandtopos.at(ligandname);
}

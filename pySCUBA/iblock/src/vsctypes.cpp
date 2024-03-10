/*
 * sidechainff.cpp
 *
 *  Created on: 2018年1月31日
 *      Author: hyliu
 */
#include "iblock/vsctypes.h"
#include "dataio/inputlines.h"
#include "dataio/datapaths.h"
#include "iblock/intrctmol.h"
#include "geometry/calculators.h"
#include "unistd.h"
using namespace NSPintrct;
std::map<std::string, int> VSCType::stericatomtypes;
std::map<char, std::string> VSCType::resnamefrom1letter;
std::vector<PackingAtomType> VSCType::packingatomtypes;
PackingEneFuncs VSCType::packingenefuncs;
PackingEneFuncs VSCType::softpackingenefuncs;
std::map<std::string, std::set<std::string>> VSCType::rotatablescatoms = { {
		"ALA", { } }, { "CYS", { "SG" } }, { "ASP", { "CG", "OD1" } }, { "GLU",
		{ "CG", "CD", "OE1" } }, { "PHE", { "CG", "CD1" } }, { "GLY", { } }, {
		"HIS", { "CG", "ND1" } }, { "ILE", { "CG1", "CD1" } }, { "LYS", { "CG",
		"CD", "CE", "NZ" } }, { "LEU", { "CG", "CD1" } }, { "MET", { "CG", "SD",
		"CE" } }, { "ASN", { "CG", "OD1" } }, { "PRO", { } }, { "GLN", { "CG",
		"CD", "OE1" } }, { "ARG", { "CG", "CD", "NE", "CZ" } }, { "SER",
		{ "OG" } }, { "THR", { "OG1" } }, { "VAL", { "CG1" } }, { "TRP", { "CG",
		"CD1" } }, { "TYR", { "CG", "CD1" } } };
const VSCType & VSCType::getVSCType(const std::string &resname,
		const std::string &filename) {
	static std::map<std::string, VSCType> vsctypes;
	static bool initialized { false };
	if (!initialized) {
#ifdef _OPENMP
#pragma omp critical(vsctype_global)
		{
			if(!initialized) {
#endif
		std::string file = filename;
		if (file.empty())
			file = std::string("sidechainff.dat");
//		file =NSPdataio::datapath()+file;
		file = NSPdataio::datafilename(file);
		vsctypes = readVSCTypes(file);
		initialized = true;
#ifdef _OPENMP
#pragma omp flush
	}
}
#endif
	}
	return vsctypes.at(resname);
}

std::map<std::string, VSCType> VSCType::readVSCTypes(
		const std::string & filename) {
	std::map<std::string, VSCType> readtypes;
	NSPdataio::InputLines inputlines;
	inputlines.init(filename, '#');
	int lidx = 0;
	std::vector<std::string> &words = inputlines[lidx++];
	int nstericatomtypes = std::stoi(words[0]);
	packingatomtypes.resize(2 * nstericatomtypes);
	for (int i = 0; i < nstericatomtypes; ++i) {
		words = inputlines[lidx++];
		int widx = 0;
		int t = std::stoi(words[widx++]);
		assert(t < 2 * nstericatomtypes);
		packingatomtypes[t].radius = std::stod(words[widx++]);
		packingatomtypes[t].hbtype = std::stoi(words[widx++]);
		packingatomtypes[t].aromatic = std::stoi(words[widx++]);
		packingatomtypes[t].nconnections = std::stoi(words[widx++]);
		for (int w = widx; w < words.size(); ++w) {
			std::string key = words[w];
			stericatomtypes.insert(std::make_pair(key, t));
		}
	}
	std::string ligandfilename=NSPdataio::datafilename("ligandatomtypes.dat");
	if(access(ligandfilename.c_str(), F_OK ) != -1){
		NSPdataio::InputLines ligandlines;
		ligandlines.init(ligandfilename,'#');
		int llnidx=0;
		std::vector<std::string> &words = ligandlines[llnidx++];
		int nlatypes=std::stoi(words[0]);
		for(int i=0;i<nlatypes;++i){
				words = ligandlines[llnidx++];
				int widx = 0;
				packingatomtypes.push_back(PackingAtomType());
				int t=packingatomtypes.size();
				packingatomtypes[t].radius = std::stod(words[widx++]);
				packingatomtypes[t].hbtype = std::stoi(words[widx++]);
				packingatomtypes[t].aromatic = std::stoi(words[widx++]);
				packingatomtypes[t].nconnections = std::stoi(words[widx++]);
				for (int w = widx; w < words.size(); ++w) {
					std::string key = words[w];
					stericatomtypes.insert(std::make_pair(key, t));
				}
		}
	}
	nstericatomtypes=stericatomtypes.size();
	packingenefuncs.setup(nstericatomtypes);
	softpackingenefuncs.setup(nstericatomtypes);
//	std::shared_ptr<EneFunc1D> hbfunc(new HarmGauEne(0.3,0.0,0.01,0.3,6));
	std::shared_ptr<EneFunc1D> hbfunc(new LJGauEne(0.29, 0.0, 0.01));
	std::shared_ptr<EneFunc1D> softhbfunc(
			new SwitchGauEne(0.3, 0.0, 0.21, 10.0, 0.14));
	for (int i = 0; i < nstericatomtypes; ++i) {
		auto & pi = packingatomtypes[i];
		bool inpolar = (pi.hbtype == 0);
		for (int j = i; j < nstericatomtypes; ++j) {
			auto & pj = packingatomtypes[j];
			if (PackingAtomType::hbond(pi.hbtype, pj.hbtype)) {
				packingenefuncs.addfunction(i, j, hbfunc);
				packingenefuncs.addfunction(j, i, hbfunc);
				softpackingenefuncs.addfunction(i, j, softhbfunc);
				softpackingenefuncs.addfunction(j, i, softhbfunc);
				continue;
			}
			bool jnpolar = (pj.hbtype == 0);
			double rmin = 0.1
					* (packingatomtypes[i].radius + packingatomtypes[j].radius);
			double emin;
			if (inpolar && jnpolar) {
				if (pi.aromatic != 0 || pj.aromatic != 0)
					emin = 1.0;
				else
					emin = 1.0;
			} else if (inpolar || jnpolar) {
				emin = 0.2;
			} else {
				emin = 0.0;
			}
			double ncon = pi.nconnections + pj.nconnections;
			if (ncon > 2)
				emin = emin * 1.5 / (ncon - 1.0);
//			std::shared_ptr<EneFunc1D> func(new HarmGauEne(rmin,emin,0.09,0.14,6));
			std::shared_ptr<EneFunc1D> func(new LJGauEne(rmin, emin, 0.14));
			packingenefuncs.addfunction(i, j, func);
			packingenefuncs.addfunction(j, i, func);
			double rcore = rmin * 0.7;
			double emax = 10.0;
			double emin_s = -1.0;
			std::shared_ptr<EneFunc1D> softfunc(
					new SwitchGauEne(rmin, emin_s, rcore, emax, 0.14));
			softpackingenefuncs.addfunction(i, j, softfunc);
			softpackingenefuncs.addfunction(j, i, softfunc);
		}
	}
	while (lidx < inputlines.size()) {
		words = inputlines[lidx++];
		std::string resname = words[0];
		readtypes.insert(std::make_pair(resname, VSCType()));
		VSCType &vsc = readtypes.at(resname);
		vsc.resname = resname;
		vsc.pdbname = words[1];
		vsc.oneletter = words[2][0];
		resnamefrom1letter.insert(std::make_pair(vsc.oneletter, vsc.resname));
		int nscatoms = std::stoi(words[3]);
		vsc.nscatoms = nscatoms;
		if (nscatoms == 0)
			continue;
		for (int i = 0; i < nscatoms; ++i) {
			int widx = 0;
			words = inputlines[lidx++];
			std::string atomname = words[widx++];
			int atype = getstericatomtype(resname, atomname);
			if (atype < 0) {
				std::cout << "Steric atomtype for " << resname << ":"
						<< atomname << " is not defined." << std::endl;
				exit(1);
			}
			double sigma = std::stod(words[widx++]);
			int hbtype = std::stoi(words[widx++]);
			assert(hbtype == packingatomtypes[atype].hbtype);
			int rotameratom = std::stoi(words[widx++]);
			int aj = std::stoi(words[widx++]);
			double b = std::stod(words[widx++]);
			double ak = std::stoi(words[widx++]);
			double a = std::stod(words[widx++]);
			double al = std::stoi(words[widx++]);
			double p = std::stod(words[widx++]);
			vsc.atomnames.push_back(atomname);
			if (rotameratom == 1)
				vsc.rotameratoms.push_back(i);
			std::vector<std::pair<int, double>> intcrd;
			intcrd.push_back(std::make_pair(aj, b));
			intcrd.push_back(std::make_pair(ak, a));
			intcrd.push_back(std::make_pair(al, p));
			vsc.internalcrds.push_back(intcrd);
		}
		words = inputlines[lidx++];
		int widx = 0;

		for (int i = 0; i < words.size() / 4; ++i) {
			int ai = std::stoi(words[widx++]);
			int aj = std::stoi(words[widx++]);
			double b = std::stod(words[widx++]);
			double b2 = std::stod(words[widx++]);
			b2 = 2 * KBT / (b2 * b2);
			vsc.newbonds.push_back(std::make_pair(ai, aj));
			vsc.b0.push_back(b);
//			vsc.kb0.push_back(b2);
			vsc.kb0.push_back(1000.0);
		}
		words = inputlines[lidx++];
		widx = 0;
//		std::cout <<resname <<std::endl;
		for (int i = 0; i < words.size() / 5; ++i) {
			std::vector<int> aijk;
			for (int m = 0; m < 3; ++m)
				aijk.push_back(std::stoi(words[widx++]));
			vsc.newangles.push_back(aijk);
			double a = std::stod(words[widx++]);
			double a2 = std::stod(words[widx++]);
			a2 = 2. * KBT * KANG_FAC / (9.0 * a2 * a2);
			vsc.a0.push_back(a);
			if (a2 > 2000.0)
				a2 = 2000.0;
			vsc.ka0.push_back(a2);
//			std::cout <<"ka "<< a2<<std::endl;
		}
		words = inputlines[lidx++];
		widx = 0;
		for (int i = 0; i < words.size() / 6; ++i) {
			std::vector<int> aijkl;
			for (int m = 0; m < 4; ++m)
				aijkl.push_back(std::stoi(words[widx++]));
			vsc.newimpdihs.push_back(aijkl);
			double p = std::stod(words[widx++]);
			double p2 = std::stod(words[widx++]);
			p2 = 2.0 * KBT * KANG_FAC / (9.0 * p2 * p2);
			vsc.imp0.push_back(p);
			vsc.kimp0.push_back(p2);
			if (p2 > 2000.0)
				p2 = 2000.0;
//			std::cout <<"kimp " <<p2 <<std::endl;
		}
		words = inputlines[lidx++];
		widx = 0;
		for (int i = 0; i < words.size() / 4; ++i) {
			std::vector<int> aijkl;
			for (int m = 0; m < 4; ++m)
				aijkl.push_back(std::stoi(words[widx++]));
			vsc.newtorsions.push_back(aijkl);
		}
	}
	return readtypes;
}
std::vector<std::vector<SCInChain>> NSPintrct::makescinchains(
		const IntrctMol &mol) {
	const std::vector<std::vector<BSInChain>> &bschains = mol.bsinchains;
	std::vector<std::vector<SCInChain>> scchains(bschains.size());
	int cidx = 0;
	for (const std::vector<BSInChain> & bsc : bschains) {
		auto & scc = scchains[cidx++];
		for (int r = 0; r < bsc.size(); ++r) {
			const IntrctBlck & blk = mol.getblck(NSPdstl::Idx2D {
					bsc[r].chainid, bsc[r].resid });
			std::vector<std::vector<int>> kaiatoms;
			const VSCType & vsc = VSCType::getVSCType(blk.resname());
			int nscatoms = vsc.nscatoms;
			assert(nscatoms == blk.natoms() - 4);
			int offset = blk.aoffset;
			for (int rsc : vsc.rotameratoms) {
				kaiatoms.push_back(std::vector<int>());
				kaiatoms.back().push_back(rsc + offset + 2);
				for (int i = 0; i < 3; ++i) {
					kaiatoms.back().push_back(
							vsc.internalcrds[rsc][i].first + offset);
				}
/*				for(int i:kaiatoms.back())
					std::cout <<mol.residuename(mol.residueofatoms()[i])
					           <<mol.atomname(i)<<":";
				std::cout<<std::endl;*/
			}
//			std::cout <<std::endl;
			scc.push_back(
					SCInChain(vsc.resname, kaiatoms, offset + 2, nscatoms));
		}
	}
	return scchains;
}
void IntrctMol::calcconformercodes(bool updateonly) const {
	if(phipsicodes.empty()) calcphipsicodes();
	if (!updateonly) {
		conformercodes.clear();
		int cidx=0;
		for (auto &sc : scinchains) {
			conformercodes.push_back(
					makeconformercodes(crds_all_, sc,phipsicodes[cidx++]));
		}
	} else {
		for (int cid = 0; cid < bsinchains.size(); ++cid) {
#pragma omp parallel for schedule(dynamic,1)
			for(int i=0;i<bsinchains[cid].size();++i){
				if(mol_(cid,i).activemod> IntrctBlck::SIDECHAIN)continue;
				if(scinchains[cid][i].kaiatoms.empty()) continue;
				ConformerCode &cc=conformercodes[cid][i];
				cc.restype=scinchains[cid][i].restype;
				cc.phipsicodes=&(phipsicodes[cid][i]);
				cc.sidechaintorsions.clear();
				cc.dtorsioncodesdx.clear();
				cc.torsioncodes.clear();
				for(auto & ijkl:scinchains[cid][i].kaiatoms){
					std::vector<NSPgeometry::XYZ> dtdx;
					double kai=NSPgeometry::torsion(crds_all_[ijkl[0]],
							crds_all_[ijkl[1]],crds_all_[ijkl[2]],crds_all_[ijkl[3]],&dtdx);
					cc.sidechaintorsions.push_back(kai);
					std::vector<DvDxi> dadx;
					dadx.push_back(std::make_pair(ijkl[0],dtdx[0]));
					dadx.push_back(std::make_pair(ijkl[1],dtdx[1]));
					dadx.push_back(std::make_pair(ijkl[2],dtdx[2]));
					dadx.push_back(std::make_pair(ijkl[3],dtdx[3]));
					std::vector<std::vector<DvDxi>> dcdx;
					cc.dtorsioncodesdx.push_back(std::vector<std::vector<DvDxi>>());
					cc.torsioncodes.push_back(
							ConformerCode::gettorsioncodes(kai,dadx,&(cc.dtorsioncodesdx.back())));
				}
			}
		}
	}
}

std::vector<double> ConformerCode::gettorsioncodes(double ang,
		const std::vector<DvDxi> &dadx, std::vector<std::vector<DvDxi>> *dcdx) {
	const std::vector<double> nang { 1.0, 2.0, 4.0 };
	std::vector<double> res;
	dcdx->clear();
	for (double n : nang) {
		double c = cos(n * ang);
		double s = sin(n * ang);
		res.push_back(c);
		dcdx->push_back(std::vector<DvDxi>());
		auto & dccdx = dcdx->back();
		for (auto &d : dadx)
			dccdx.push_back(-(s * n) * d);
		res.push_back(s);
		dcdx->push_back(std::vector<DvDxi>());
		auto & dcsdx = dcdx->back();
		for (auto &d : dadx)
			dcsdx.push_back((c * n) * d);
	}
	return res;
}
std::vector<ConformerCode> NSPintrct::makeconformercodes(
		const std::vector<NSPgeometry::XYZ> &crds,
		const std::vector<SCInChain> &scinchains,
		std::vector<PhiPsiCodes> &phipsicodes) {
	int nscs = scinchains.size();
	std::vector<ConformerCode> result(nscs);
#pragma omp parallel for schedule(dynamic,1)
	for (int i = 0; i < nscs; ++i) {
		if (scinchains[i].kaiatoms.empty())
			continue;
		ConformerCode &cc = result[i];
		cc.restype = scinchains[i].restype;
		cc.phipsicodes = &(phipsicodes[i]);
		for (auto & ijkl : scinchains[i].kaiatoms) {
			std::vector<NSPgeometry::XYZ> dtdx;
			double kai = NSPgeometry::torsion(crds[ijkl[0]], crds[ijkl[1]],
					crds[ijkl[2]], crds[ijkl[3]], &dtdx);
			cc.sidechaintorsions.push_back(kai);
			std::vector<DvDxi> dadx;
			dadx.push_back(std::make_pair(ijkl[0], dtdx[0]));
			dadx.push_back(std::make_pair(ijkl[1], dtdx[1]));
			dadx.push_back(std::make_pair(ijkl[2], dtdx[2]));
			dadx.push_back(std::make_pair(ijkl[3], dtdx[3]));
			std::vector<std::vector<DvDxi>> dcdx;
			cc.dtorsioncodesdx.push_back(std::vector<std::vector<DvDxi>>());
			cc.torsioncodes.push_back(
					ConformerCode::gettorsioncodes(kai, dadx,
							&(cc.dtorsioncodesdx.back())));
		}
	}
	return result;
}

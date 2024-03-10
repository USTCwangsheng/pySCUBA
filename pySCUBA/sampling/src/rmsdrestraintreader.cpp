/*
 * rmsdrestraintreader.cpp
 *
 *  Created on: 2020年9月30日
 *      Author: hyiu
 */
#include "sampling/rmsdrestraintreader.h"
#include "sampling/molsdrun.h"
#include "dataio/inputlines.h"
#include "iblock/blkselector.h"
#include "proteinrep/aaconformer.h"
using namespace NSPsampling;
enum {
	ALL_ALL, SEL_SEL
};
enum {
	MCATOMS, CAATOMS, ALLATOMS
};
const std::array<std::string, 2> blkset { "all_all", "sel_sel" };
const std::array<std::string, 3> atomset { "mcatoms", "caatoms", "allatoms" };
/*void RMSDRestraintReader::addstructrestraints(MolSDRun &sdrun) {
	if (!sdrun.initpara().rmsdrestraintfile.empty()) {
		readterms(sdrun.initpara().rmsdrestraintfile);
		for (auto & t : terms_)
			addrestraint(sdrun, t);
	}
}
void RMSDRestraintReader::addrestraint(MolSDRun &sdrun, RMSDTerm &term) {
	NSPproteinrep::AAConformersInModel refmodel;
	refmodel.readpdbfile(term.pdbref);
	auto &refconformers = refmodel.conformers;
	std::vector<NSPgeometry::XYZ> refcrds = sdrun.imol()->recollectcrds_all();
	const NSPintrct::IntrctMol & im = *(sdrun.imol());
	std::vector<double> weight(refcrds.size(), 0.0);
	if (term.refblks.empty()) {
		for (int cid = 0; cid < im.nchains(); ++cid) {
			for (int bid = 0; bid < im.nresidues(cid); ++bid) {
				term.refblks.push_back(std::make_pair(cid, bid));
				term.restrainedblks.push_back(std::make_pair(cid, bid));
			}
		}
	}
	for (int i = 0; i < term.refblks.size(); ++i) {
		NSPdstl::Idx2D blkidx(term.restrainedblks[i].first,
				term.restrainedblks[i].second);
		NSPproteinrep::AAConformer & refaa =
				refconformers[term.refblks[i].first][term.refblks[i].second];
		if (term.restrainedatoms == CAATOMS) {
			int caidx = im.getblck(blkidx).aoffset + 1;
			refcrds[caidx] = refaa.getglobalcrd().at("CA") * A2NM;
			weight[caidx] = 1.0;
		} else if (term.restrainedatoms == MCATOMS) {
			int nidx = im.getblck(blkidx).aoffset;
			refcrds[nidx] = refaa.getglobalcrd().at("N") * A2NM;
			int caidx = im.getblck(blkidx).aoffset + 1;
			refcrds[caidx] = refaa.getglobalcrd().at("CA") * A2NM;
			int cidx = im.getblck(blkidx).aoffset + im.getblck(blkidx).natoms()
					- 2;
			refcrds[cidx] = refaa.getglobalcrd().at("C") * A2NM;
			int oidx = im.getblck(blkidx).aoffset + 1
					+ im.getblck(blkidx).natoms() - 1;
			refcrds[oidx] = refaa.getglobalcrd().at("O") * A2NM;
			weight[nidx] = 1.0;
			weight[caidx] = 1.0;
			weight[cidx] = 1.0;
			weight[oidx] = 1.0;
		} else if (term.restrainedatoms == ALLATOMS) {
			for (int a = 0; a < im.getblck(blkidx).natoms(); ++a) {
				std::string atomname = im.getblck(blkidx).topo->atomnames.at(a);
				refcrds[im.getblck(blkidx).aoffset + a] =
						refaa.getglobalcrd().at(atomname) * A2NM;
				weight[im.getblck(blkidx).aoffset + a] = 1.0;
			}
		}
	}
	auto newrestr = std::shared_ptr<ResTraint>(
			new StructRestraint(refcrds, weight, term.kres / (A2NM * A2NM),
					term.rmsdmin * A2NM, term.rmsdswitch * A2NM,
					StructRestraint::TOTALMODE));
	sdrun.addrestraint(newrestr);
}*/
void RMSDRestraintReader::addrestraints(NSPintrct::IntrctMol &imol,
		std::vector<std::shared_ptr<ResTraint>> &results, int& switchkey) {
	for (auto & term : terms_) {
		NSPproteinrep::AAConformersInModel refmodel;
		refmodel.readpdbfile(term.pdbref);
		auto &refconformers = refmodel.conformers;
		std::vector<NSPgeometry::XYZ> refcrds = imol.recollectcrds_all();
		const NSPintrct::IntrctMol & im = imol;
		std::vector<double> weight(refcrds.size(), 0.0);
		if (term.refblks.empty()) {
			for (int cid = 0; cid < im.nchains(); ++cid) {
				for (int bid = 0; bid < im.nresidues(cid); ++bid) {
					term.refblks.push_back(std::make_pair(cid, bid));
					term.restrainedblks.push_back(std::make_pair(cid, bid));
				}
			}
		}
		for (int i = 0; i < term.refblks.size(); ++i) {
			NSPdstl::Idx2D blkidx(term.restrainedblks[i].first,
					term.restrainedblks[i].second);
			NSPproteinrep::AAConformer & refaa =
					refconformers[term.refblks[i].first][term.refblks[i].second];
			if (term.restrainedatoms == CAATOMS) {
				int caidx = im.getblck(blkidx).aoffset + 1;
				refcrds[caidx] = refaa.getglobalcrd().at("CA") * A2NM;
				weight[caidx] = 1.0;
			} else if (term.restrainedatoms == MCATOMS) {
				int nidx = im.getblck(blkidx).aoffset;
				refcrds[nidx] = refaa.getglobalcrd().at("N") * A2NM;
				int caidx = im.getblck(blkidx).aoffset + 1;
				refcrds[caidx] = refaa.getglobalcrd().at("CA") * A2NM;
				int cidx = im.getblck(blkidx).aoffset
						+ im.getblck(blkidx).natoms() - 2;
				refcrds[cidx] = refaa.getglobalcrd().at("C") * A2NM;
				int oidx = im.getblck(blkidx).aoffset + 1
						+ im.getblck(blkidx).natoms() - 1;
				refcrds[oidx] = refaa.getglobalcrd().at("O") * A2NM;
				weight[nidx] = 1.0;
				weight[caidx] = 1.0;
				weight[cidx] = 1.0;
				weight[oidx] = 1.0;
			} else if (term.restrainedatoms == ALLATOMS) {
				for (int a = 0; a < im.getblck(blkidx).natoms(); ++a) {
					std::string atomname =
							im.getblck(blkidx).topo->atomnames.at(a);
					refcrds[im.getblck(blkidx).aoffset + a] =
							refaa.getglobalcrd().at(atomname) * A2NM;
					weight[im.getblck(blkidx).aoffset + a] = 1.0;
				}
			}
		}
		auto newrestr = std::shared_ptr<ResTraint>(
				new StructRestraint(refcrds, weight, term.kres / (A2NM * A2NM),
						term.rmsdmin * A2NM, term.rmsdswitch * A2NM,
						StructRestraint::TOTALMODE));
		newrestr->switchable() = switchkey;
		newrestr->restriantname() = "RMSD";
		results.push_back(newrestr);
	}
}
void RMSDRestraintReader::readterms(NSPdataio::InputLines &input, int &lidx) {
	while (lidx < input.size()) {
		auto & words = input[lidx++];
		if (words[0].substr(0, 3) == "END")
			break;
		if (words[0] != "RMSDterms")
			continue;
		double kres = std::stod(words[1]);
		double rmsdmin = std::stod(words[2]);
		double rmsdswitch = std::stod(words[3]);
		int rblkset { -1 };
		int idx = 0;
		for (auto & m : blkset) {
			if (words[4] == m)
				rblkset = idx;
			idx++;
		}
		assert(rblkset >= 0);
		int aset { -1 };
		idx = 0;
		for (auto & m : atomset) {
			if (words[5] == m)
				aset = idx;
			idx++;
		}
		assert(aset >= 0);
		std::string refpdbfile = words[6];
		std::vector<std::pair<int, int>> refblks;
		std::vector<std::pair<int, int>> restrainedblks;
		if (rblkset == SEL_SEL) {
			auto & w1 = input[lidx++];
			assert(w1[0] == "refresidues");
			std::string sele1;
			for (int m = 1; m < w1.size(); ++m)
				sele1 = sele1 + " " + w1[m];
			NSPintrct::BlkSelector sel1(sele1);
			for (auto & bc : sel1.selectedblks) {
				for (auto &b : bc.second) {
					refblks.push_back(std::make_pair(bc.first, b));
				}
			}
			auto & w2 = input[lidx++];
			assert(w2[0] == "restrainedresidues");
			std::string sele2;
			for (int m = 1; m < w2.size(); ++m)
				sele2 = sele2 + " " + w2[m];
			NSPintrct::BlkSelector sel2(sele2);
			for (auto & bc : sel2.selectedblks) {
				for (auto &b : bc.second) {
					restrainedblks.push_back(std::make_pair(bc.first, b));
				}
			}
			if (refblks.size() != restrainedblks.size()) {
				std::cout << "refresidues must be equal to restrainedresidues" << std::endl;
				exit(1);
			}
		}
		terms_.push_back(RMSDRestraintReader::RMSDTerm());
		auto &term = terms_.back();
		term.kres = kres;
		term.rmsdmin = rmsdmin;
		term.rmsdswitch = rmsdswitch;
		term.refblks = refblks;
		term.restrainedblks = restrainedblks;
		term.restrainedatoms = aset;
		term.pdbref = refpdbfile;
	}
}


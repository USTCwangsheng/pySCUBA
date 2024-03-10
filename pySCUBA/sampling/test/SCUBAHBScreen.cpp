/*
 * SCUBAHBScreen.cpp
 *
 *  Created on: 2021年6月3日
 *      Author: yxchen
 */


#include "iblock/molmodeler.h"
#include "iblock/intrctmol.h"
#include "iblock/intrctblck.h"
#include "sampling/screen.h"
#include <fstream>
#include <sstream>
#include <dirent.h>
#include <map>
#include "proteinrep/aaconformer.h"

using namespace std;
using namespace NSPintrct;
using namespace Screen;
using namespace NSPproteinrep;

// screen/rank scaffold(s) by HB conditions.

int main(int argc, char** argv)
{
	if (argc != 2)
	{
		cout << "Usage: SCUBAHBScreen parfile" << endl;
		exit(1);
	}
	HBPar hbp(argv[1]);
	map<double, vector<string>> nhb_sca; // not exposed & no-HB_percentage for each scaffold.
	map<string, vector<Screen::Atom>> sca_patms;
	map<string, vector<Screen::Atom>> sca_nhbatms;
	ifstream ifsd("expose_detailed.txt");
	if (ifsd.good()) system("rm expose_detailed.txt");
	ofstream ofsd("expose_detailed.txt", ios::app);
	for (auto sca : hbp.scaffolds)
	{
		ofsd << sca << ":" << endl;
		int nexp_patm = 0; // non-exposed polar atm.
		int nexp_nhbatm = 0; // non-exposed polar atm without forming hb.
		vector<Screen::Atom> patms;
		vector<Screen::Atom> nhbatms;
		MolModeler modeler;
		AAConformersInModel model;
		model.readpdbfile(sca);
		if (hbp.objectpair[0] == "MC")
		{
			auto &aaconfs = model.conformers;
			for (auto &c : aaconfs)
				for (auto &r : c)
				{
					bool inregion = false;
					if (hbp.regions.size() == 1 && hbp.regions[0].cid == "-1"
							&& hbp.regions[0].rbeg == -1 && hbp.regions[0].rend == -1)
						inregion = true;
					else
					{
						for (auto reg : hbp.regions)
						{
							if (r.chainid_or == reg.cid[0] && r.residueid_or >= reg.rbeg && r.residueid_or <= reg.rend)
								inregion = true;
							if (inregion) break;
						}
					}
					if (inregion)
					{
						if (!expose(r.globalcrd["N"], string {r.chainid_or}, r.residueid_or, aaconfs, hbp.expose_th, hbp.objectpair[1]))
						{
							nexp_patm++;
							Screen::Atom a(string {r.chainid_or}, r.residueid_or, "N");
							patms.push_back(a);
							if (!formhb(r.globalcrd["N"], string {r.chainid_or}, r.residueid_or, aaconfs, hbp.hb_th, hbp.objectpair[1]))
							{
								nexp_nhbatm++;
								nhbatms.push_back(a);
							}
						}
						if (!expose(r.globalcrd["O"], string {r.chainid_or}, r.residueid_or, aaconfs, hbp.expose_th, hbp.objectpair[1]))
						{
							nexp_patm++;
							Screen::Atom a(string {r.chainid_or}, r.residueid_or, "O");
							patms.push_back(a);
							if (!formhb(r.globalcrd["O"], string {r.chainid_or}, r.residueid_or, aaconfs, hbp.hb_th, hbp.objectpair[1]))
							{
								nexp_nhbatm++;
								nhbatms.push_back(a);
							}
						}
					}
				}
		}
		else if (hbp.objectpair[0] == "all")
		{
			auto &aaconfs = model.conformers;
			for (auto &c : aaconfs)
				for (auto &r : c)
				{
					bool inregion = false;
					if (hbp.regions.size() == 1 && hbp.regions[0].cid == "-1"
							&& hbp.regions[0].rbeg == -1 && hbp.regions[0].rend == -1)
						inregion = true;
					else
					{
						for (auto reg : hbp.regions)
						{
							if (r.chainid_or == reg.cid[0] && r.residueid_or >= reg.rbeg && r.residueid_or <= reg.rend)
								inregion = true;
							if (inregion) break;
						}
					}
					if (inregion)
					{
						for (auto it = r.globalcrd.begin(); it != r.globalcrd.end(); it++)
						{
							if (it->first[0] == 'N' || it->first[0] == 'O')
							{
								if (!expose(it->second, string {r.chainid_or}, r.residueid_or, aaconfs, hbp.expose_th, hbp.objectpair[1]))
								{
									nexp_patm++;
									Screen::Atom a(string {r.chainid_or}, r.residueid_or, it->first);
									patms.push_back(a);
									if (!formhb(it->second, string {r.chainid_or}, r.residueid_or, aaconfs, hbp.hb_th, hbp.objectpair[1]))
									{
										nexp_nhbatm++;
										nhbatms.push_back(a);
									}
								}
							}
						}
					}
				}
		}
		else if (hbp.objectpair[0] == "SC")
		{
			auto &aaconfs = model.conformers;
			for (auto &c : aaconfs)
				for (auto &r : c)
				{
					bool inregion = false;
					if (hbp.regions.size() == 1 && hbp.regions[0].cid == "-1"
							&& hbp.regions[0].rbeg == -1 && hbp.regions[0].rend == -1)
						inregion = true;
					else
					{
						for (auto reg : hbp.regions)
						{
							if (r.chainid_or == reg.cid[0] && r.residueid_or >= reg.rbeg && r.residueid_or <= reg.rend)
								inregion = true;
							if (inregion) break;
						}
					}
					if (inregion)
					{
						for (auto scatm : r.sidechainatoms[r.residuename])
							if (scatm[0] == 'N' || scatm[0] == 'O')
							{
								if (!expose(r.globalcrd[scatm], string {r.chainid_or}, r.residueid_or, aaconfs, hbp.expose_th, hbp.objectpair[1]))
								{
									nexp_patm++;
									Screen::Atom a(string {r.chainid_or}, r.residueid_or, scatm);
									patms.push_back(a);
									if (!formhb(r.globalcrd[scatm], string {r.chainid_or}, r.residueid_or, aaconfs, hbp.hb_th, hbp.objectpair[1]))
									{
										nexp_nhbatm++;
										nhbatms.push_back(a);
									}
								}
							}
					}
				}
		}
		double nhb = double(nexp_nhbatm)/nexp_patm;
		if (nhb_sca.count(nhb) > 0)
			nhb_sca[nhb].push_back(sca);
		else
			nhb_sca.insert(make_pair(nhb, vector<string> {sca}));
		sca_patms.insert(make_pair(sca, patms));
		sca_nhbatms.insert(make_pair(sca, nhbatms));
	}
	ofsd.close();

	// output
	ofstream ofs("rankbynohb.txt");
	for (auto iter = nhb_sca.begin(); iter != nhb_sca.end(); iter++)
		for (auto scafile : iter->second)
			ofs <<  iter->first << " " << scafile << endl;
	ofs.close();
	ofstream ofse("atoms_notexposed.txt");
	for (auto iter = sca_patms.begin(); iter != sca_patms.end(); iter++)
	{
		ofse << iter->first << ":" << endl;
		for (auto atm : iter->second)
			ofse << " " << atm.cid << " " << atm.rid << " " << atm.aname << endl;
	}
	ofse.close();
	ofstream ofsh("notexposed_nohb_atoms.txt");
	for (auto iter = sca_nhbatms.begin(); iter != sca_nhbatms.end(); iter++)
	{
		ofsh << iter->first << ":" << endl;
		for (auto atm : iter->second)
			ofsh << " " << atm.cid << " " << atm.rid << " " << atm.aname << endl;
	}
	ofsh.close();

}

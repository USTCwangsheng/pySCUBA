/*
 * testisolatedpocket.cpp
 *
 *  Created on: 2020年7月13日
 *      Author: hyiu
 */
#include "proteinrep/isolatedpocket.h"
#include "proteinrep/rotamerutil.h"
#include <fstream>
using namespace NSPproteinrep;

int main() {
	IsolatedPocket pocket;
	std::cout << "Read components from an existing PDB file?[Y/N]";
	std::string ans;
	getline(std::cin, ans);
	if (ans[0] == 'Y' || ans[0] == 'y') {
		std::cout << "Specify the input PDB filename: ";
		std::string pdbfilename;
		getline(std::cin, pdbfilename);
		PdbReader pdb;
		try {
			pdb.readpdb(pdbfilename);
		} catch (std::exception &e) {
			std::cout << "Error occurred reading pdbfile " << pdbfilename
					<< std::endl;
			exit(0);
		}
		std::cout << "Specify residues from the pdbfile to include in the pocket below."
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
				std::cin.getline(str, 50);
				itext = std::string(str);
				if (itext.empty())
					break;
				cmpntids.push_back(IsolatedPocket::ComponentID(itext));
			} catch (std::exception &e) {
				;
			}
		}
		pocket.addcomponents(pdb, cmpntids);
	}

	while (true) {
		std::cout << "Add another extra residue as components?[Y/N]";
		getline(std::cin, ans);
		if (ans[0] != 'Y' && ans[0] != 'y')
			break;
		char newchainid;
		std::cout << "Specify its chain ID: ";
		std::cin >> newchainid;
		getline(std::cin, ans);
		std::cout << "Specify its residue type(capital three-letter) : ";
		std::string resname;
		getline(std::cin, resname);
		AAConformer newresidue = make_aaconformer(randomrotamer(resname));
		if (pocket.empty()) {
			pocket.addcomponent(newresidue, newchainid,
					std::vector<IsolatedPocket::AtomID>(),
					std::vector<std::string>(), std::vector<double>());
			continue;
		}
		std::vector<IsolatedPocket::AtomID> anchors1(3);
		for (int i = 0; i < 3; i++) {
			std::cout << "Pocket anchor atom " << i + 1 << " : ";
			std::cin >> anchors1[i].cmpntID.chainid
					>> anchors1[i].cmpntID.residueid
					>> anchors1[i].cmpntID.resname >> anchors1[i].atmnm;
			getline(std::cin, ans);
		}
		std::vector<std::string> anchors2(3);
		std::vector<double> intcrds(6);
		std::cout << "Residue anchor atom 1: ";
		getline(std::cin, anchors2[0]);
		std::cout << "Distance to pocket anchor atom 3: ";
		std::cin >> intcrds[2];
		getline(std::cin, ans);
		std::cout << "Angle with pocket anchor atoms 3 and 2: ";
		std::cin >> intcrds[1];
		getline(std::cin, ans);
		std::cout << "Torsion with  pocket anchor atoms 3, 2 and 1: ";
		std::cin >> intcrds[0];
		getline(std::cin, ans);
		std::cout << "Residue anchor atom 2: ";
		getline(std::cin, anchors2[1]);
		std::cout << "Angle with atom 1 and pocket anchor atom 3: ";
		std::cin >> intcrds[4];
		getline(std::cin, ans);
		std::cout << "Torsion with atom 1, pocket anchor atoms 3 and 2 ";
		std::cin >> intcrds[3];
		getline(std::cin, ans);
		std::cout << "Residue anchor atom 3: ";
		getline(std::cin, anchors2[2]);
		std::cout << "Torsion with atoms 2, 1 and pocket anchor atoms 3 and 2 ";
		std::cin >> intcrds[5];
		getline(std::cin, ans);
		double toradius = 3.14159265 / 180.0;
		for (auto & v : intcrds)
			v *= toradius;
		intcrds[2] = intcrds[2] / toradius;
		pocket.addcomponent(newresidue, newchainid, anchors1, anchors2,
				intcrds);
	}
	std::cout <<"Specify output pdb file: ";
	std::string pocketpdb;
	getline(std::cin, pocketpdb);
	std::ofstream ofs(pocketpdb);
	pocket.writepdb(ofs);
}

#include "iblock/molmodeler.h"
#include "iblock/intrctmol.h"
#include "iblock/intrctblck.h"
#include "designseq/StructureInfo.h"
#include "geometry/calculators.h"
#include "designseq/ProteinRep.h"
#include "dstl/randomengine.h"
#include "proteinrep/residuestate.h"
#include "proteinrep/aminoacidseq.h"

#include <dataio/splitstring.h>
#include <fstream>
#include <cassert>
#include <string>
#include <iostream> 
/*��ȡpdb�����ṹ,��ȡloop��*/
using namespace std;
using namespace NSPintrct;
using namespace NSPdesignseq;
using namespace NSPproteinrep;
using namespace NSPgeometry;

int main(int argc, char** argv) {
	if (argc != 2) {
		std::cout << "usage: readss pdb" << std::endl;
		exit(0);
	}
	std::string pdbfile(argv[1]);
	std::ifstream is(pdbfile);
	auto startrstates = residuestates(is);
	int chainidx = 0;
	std::string startss;
	for (auto& c : startrstates) {
		std::string seq;
		int seqidx = 1;
		for (auto& r : c) {
			seq.push_back(AminoAcidSeq::name2code(r.sidechainstate.resiudetype));
			startss.push_back(r.backbonestate.SSState);
		}
		//std::cout << seq << std::endl;
		std::cout<< startss<<std::endl;
	}
	std::vector <string> cssn;
	int last=0;
	char buffer[100];
	for (int i = 2;i < startss.size()-2;i++) {
		if (startss[i] == 'C') {
			if (startss[i - 1] != 'C' || startss[i + 1] != 'C') {
				
				if (cssn.size() % 4 == 0) {
					cssn.push_back("A ");
				}
				if (startss[i - 1] != 'C') {
					if(startss[i - 1] == 'H' && startss[i - 5] == 'H'){ 
						cssn.push_back(to_string(i + 1 - 1) + " "); 
					}
					else if (startss[i - 1] == 'E') {
						cssn.push_back(to_string(i + 1 - 1) + " ");
					}
				}
				if (startss[i + 1] != 'C') {
					if (startss[i + 1] == 'H' && startss[i + 5] == 'H') {
						//std::cout << startss[i + 5] << std::endl;
						cssn.push_back(to_string(i + 1 + 1) + " ");
					}
					else if (startss[i + 1] == 'E'){
						cssn.push_back(to_string(i + 1 + 1) + " ");
					}
				}
				if (cssn.size() % 4 == 3) {
					istringstream ss(cssn[cssn.size()-2]);
					int t;
					ss >> t;
					std::cout << t << std::endl;
					if(i - t + 1<4){ cssn.pop_back(); cssn.pop_back();
					}
					else if(i - t + 1>7){ cssn.push_back(to_string(i - t) + " "); }
					else cssn.push_back(to_string(i - t + 1)+" ");
					
				}
				
				
				std::cout << i << std::endl;
			}
		}
	
	}
	for (auto& n : cssn) { std::cout << n << std::endl; }
	std::string filename = "loopidx.txt";
	std::ofstream ofs1(filename);
	for (auto c : cssn) ofs1 << c;
	ofs1 << std::endl;

	std::string filename2 = "ss.txt";
	std::ofstream ofs2(filename2);
	for (auto n : startss) ofs2 << n;
	ofs2 << std::endl;
}
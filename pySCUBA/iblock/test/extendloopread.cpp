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
	if (argc != 3) {
		std::cout << "usage: readss pdb intendedssfile" << std::endl;
		exit(0);
	}
	std::string pdbfile(argv[1]);
	std::string intendssfile(argv[2]);
	std::ifstream is(pdbfile);
	std::ifstream ssfileline(intendssfile);
	std::string intenedss;
	getline(ssfileline, intenedss);
	std::cout << intenedss << std::endl;
	auto startrstates = residuestates(is);
	int chainidx = 0;
	std::string startss, readss;
	for (auto& c : startrstates) {
		std::string seq;
		int seqidx = 1;
		for (auto& r : c) {
			seq.push_back(AminoAcidSeq::name2code(r.sidechainstate.resiudetype));
			readss.push_back(r.backbonestate.SSState);
		}
		//std::cout << seq << std::endl;
		std::cout<< readss <<std::endl;
	}
	if (intenedss.length() == readss.length()) {
		for (int i = 0;i < intenedss.size();i++) {
			if (readss[i] == 'E') { startss.push_back(readss[i]); }
			else startss.push_back(intenedss[i]);
		}
	}
	else { std::cout << "intenedss.length() != readss.length()" << std::endl; exit(0); }
	assert(startss.length() == readss.length());
	std::cout << startss << std::endl;
	std::vector <string> cssn;
	int last=0;
	char buffer[100];
	for (int i = 2;i < startss.size()-2;i++) {
		if (startss[i] == 'C') {
			if (startss[i - 1] != 'C' || startss[i + 1] != 'C') {
				
				if (cssn.size() % 4 == 0) {
					cssn.push_back("A ");
				}
				if(cssn.size() % 4 != 3) {
					if (startss[i - 1] != 'C' && startss[i + 1] != 'C') {
						if (startss[i - 1] == 'H' && startss[i - 2] == 'H' && startss[i - 3] == 'H' && startss[i - 4] == 'H' && startss[i - 5] == 'H'){ cssn.push_back(to_string(i + 1 + 1) + " "); }
						else if (startss[i + 1] == 'H' && startss[i + 2] == 'H' && startss[i + 3] == 'H' && startss[i + 4] == 'H' && startss[i + 5] == 'H') { cssn.push_back(to_string(i + 1 - 1) + " "); }
						else {
							if (cssn.size() % 4 == 1)
							{
								cssn.push_back(to_string(i + 1 - 1) + " "); cssn.push_back(to_string(i + 1 + 1) + " ");
							}
							else { cssn.push_back(to_string(i + 1 + 1) + " "); }
						}
					}
					if (startss[i - 1] != 'C') {
					if(startss[i - 1] == 'H'&& startss[i - 2] == 'H'&& startss[i - 3] == 'H'&& startss[i - 4] == 'H' && startss[i - 5] == 'H'){
						cssn.push_back(to_string(i + 1 - 1) + " "); 
					}
					else if (startss[i - 1] == 'E') {
						cssn.push_back(to_string(i + 1 - 1) + " ");
					}
				}
				else if (startss[i + 1] != 'C') {
					if (startss[i + 1] == 'H' && startss[i + 2] == 'H'&& startss[i + 3] == 'H'&& startss[i + 4] == 'H' && startss[i + 5] == 'H') {
						//std::cout << startss[i + 5] << std::endl;
						cssn.push_back(to_string(i + 1 + 1) + " ");
					}
					else if (startss[i + 1] == 'E'){
						cssn.push_back(to_string(i + 1 + 1) + " ");
					}
				}
				}
				if (cssn.size() % 4 == 3) {
					istringstream ss2(cssn[cssn.size()-2]);//�����ڶ���Ԫ��,����תΪint
					int t2;
					ss2 >> t2;
					//std::cout << t << std::endl;
					if (i - t2 + 1 == 1) {
						if (startss[i + 1] == 'E') {
							cssn[cssn.size() - 2] = to_string(t2 - 2) + " ";//��ǰ��H-2
							cssn[cssn.size() - 1] = to_string(i + 2 + 1) + " ";//�����E+1
						}
						else {
							cssn[cssn.size() - 2] = to_string(t2 - 1) + " ";//��ǰ��E-1
							cssn[cssn.size() - 1] = to_string(i + 2 + 2) + " ";//�����H+2
						}
						cssn.push_back(to_string(4) + " ");
					}
					else if (i - t2 + 1 == 2) {
						if (startss[i + 1] == 'E') {
							cssn[cssn.size() - 2] = to_string(t2 - 2) + " ";//��ǰ��H-2
						}
						else {
							cssn[cssn.size() - 1] = to_string(i + 2 + 2) + " ";//�����H+2
						}
						cssn.push_back(to_string(4) + " ");
					}
					else if (i - t2 + 1 == 3) {
						if (startss[i + 1] == 'E') {
							cssn[cssn.size() - 2] = to_string(t2 - 1) + " ";//��ǰ��H-1
						}
						else {
							cssn[cssn.size() - 1] = to_string(i + 2 + 1) + " ";//�����H+1
						}
						cssn.push_back(to_string(4) + " ");
					}
					//else if(i - t2 + 1>6){ cssn.push_back(to_string(i - t2) + " "); }
					else if (i - t2 + 1 ==5) { 
						if (startss[i + 1] == 'E') {
							cssn[cssn.size() - 2] = to_string(t2 + 0) + " ";//��ǰ��H+0
							cssn.push_back(to_string(5) + " ");
						}
						else {
							cssn[cssn.size() - 1] = to_string(i + 2 - 1) + " ";//�����H-1
							cssn.push_back(to_string(4) + " ");
						}
					}
					else if (i - t2 + 1 == 6) {
						if (startss[i + 1] == 'E') {
							cssn[cssn.size() - 2] = to_string(t2 + 1) + " ";//��ǰ��H+1
							cssn.push_back(to_string(5) + " ");
						}
						else {
							cssn[cssn.size() - 1] = to_string(i + 2 - 2) + " ";//�����H-2
							cssn.push_back(to_string(4) + " ");
						}
						
					}
					else if (i - t2 + 1 > 6) {
						if (startss[i + 1] == 'E') {
							cssn[cssn.size() - 2] = to_string(t2 + 0) + " ";//��ǰ��H+0
							cssn.push_back(to_string(i - t2 + 0) + " ");//���ȼ�1
						}
						else {
							cssn[cssn.size() - 1] = to_string(i + 2 - 2) + " ";//�����H-2
							cssn.push_back(to_string(i - t2 - 1) + " ");//���ȼ�2
						}
					}
					else cssn.push_back(to_string(i - t2 + 1)+" ");//4��ԭ��
					
				}
				
				
				std::cout << i << std::endl;
			}
		}
		if (cssn.size() % 4 != 2 && cssn.size() % 4 != 3 && startss[i] == 'E' && startss[i + 1] == 'H' && startss[i + 2] == 'H' && startss[i + 3] == 'H' && startss[i + 4] == 'H' && startss[i + 5] == 'H') {
			cssn.push_back("A ");cssn.push_back(to_string(i + 1 - 1) + " "); cssn.push_back(to_string(i + 1 + 4) + " ");cssn.push_back(to_string(4) + " ");
		}
		if (cssn.size() % 4 != 2 && cssn.size() % 4 != 3 && startss[i] == 'H' && startss[i + 1] == 'E' && startss[i - 1] == 'H' && startss[i - 2] == 'H' && startss[i - 3] == 'H' && startss[i - 4] == 'H' && startss[i - 5] == 'H') {
			cssn.push_back("A ");cssn.push_back(to_string(i + 1 - 4) + " "); cssn.push_back(to_string(i + 1 + 1) + " ");cssn.push_back(to_string(4) + " ");
		}
	}
	for (auto& n : cssn) { std::cout << n << std::endl; }
	std::string filename = "loopidx.txt";
	int loopnumber = cssn.size() / 4;
	std::ofstream ofs1(filename);
	for (auto c : cssn) ofs1 << c;
	ofs1 << std::endl;
	ofs1 << loopnumber << std::endl;
	ofs1.close();

	std::string filename2 = "ss.txt";
	std::ofstream ofs2(filename2);
	for (auto n : startss) ofs2 << n;
	ofs2 << std::endl;
	ofs2.close();
}
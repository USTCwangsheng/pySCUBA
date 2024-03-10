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

using namespace std;
using namespace NSPintrct;
using namespace NSPdesignseq;
using namespace NSPproteinrep;
using namespace NSPgeometry;

int main(int argc, char** argv) {
	if (argc != 3) {
		std::cout << "usage: comparess startpdb topn" << std::endl;
		exit(0);
	}
	std::string pdbfile(argv[1]);
	std::ifstream is(pdbfile);
	std::string topnstring(argv[2]);
	int topn=std::stoi(topnstring);
	auto startrstates = residuestates(is);
	int chainidx = 0;
	std::string startss;

	int start_H_num = 0;
	int start_E_num = 0;
	int start_C_num = 0;
	for (auto& c : startrstates) {
		std::string seq;

		//std::cout << "Chain " << chainidx++ << ":" << std::endl;
		//std::cout << "seq\tresidue\tSS\tSAI\tPHI\tPSI" << std::endl;
		int seqidx = 1;
		for (auto& r : c) {
			seq.push_back(AminoAcidSeq::name2code(r.sidechainstate.resiudetype));
			startss.push_back(r.backbonestate.SSState);
			//	std::cout << seqidx++ << "\t"
			//		<< r.sidechainstate.resiudetype << "\t"
			//		<< r.backbonestate.SSState << "\t"
			//		<< r.backbonestate.sai << "\t"
			//		<< r.backbonestate.phi << "\t"
			//		<< r.backbonestate.psi << std::endl;
			//}
		}
		//std::cout << seq << std::endl;
		std::cout << startss << std::endl;
        
		std::vector<std::string> words;
		NSPutils::split(startss,words,"CE");
		for (auto &word : words){
            std:cout << word << endl;
			start_H_num = start_H_num + word.size();
		}
        std::cout << "num of startH" << start_H_num << std::endl;
		
        words.clear();
		NSPutils::split(startss,words,"CH");
		for (auto &word : words){
            std::cout << word << endl;
			start_E_num = start_E_num + word.size();
		}
        std::cout << "num of startE" << start_E_num << std::endl;

        words.clear();
		NSPutils::split(startss,words,"EH");
		for (auto &word : words){
            std::cout << word << endl;
			start_C_num = start_C_num + word.size();
		}
        std::cout << "num of startC" << start_C_num << std::endl;
	}
	double av_H_num;
	double av_E_num;
	double av_C_num;
	int tot_H_num=0;
	int tot_E_num=0;
	int tot_C_num=0;

	for (int i = 0;i <= topn;++i) {
		std::string topn_name;
		topn_name= "enetop_" + std::to_string(i) + ".pdb";
		std::ifstream topnfile(topn_name);
		if (topnfile.good()==false) { 
			std::cout << "no enetop_n file" << std::endl;
			exit(1); }
		auto topnrstates = residuestates(topnfile);
		std::string topnss;

		for (auto& c : topnrstates) {
			std::string topnseq;

			int seqidx = 1;
			for (auto& r : c) {
				topnseq.push_back(AminoAcidSeq::name2code(r.sidechainstate.resiudetype));
				topnss.push_back(r.backbonestate.SSState);

			}
			//std::cout << seq << std::endl;
			//std::cout << topnss << std::endl;

			int H_num = 0;
			std::vector<std::string> words;
		    NSPutils::split(topnss,words,"CE");
			for (auto &word : words){
				H_num = H_num + word.size();
			}
			tot_H_num = tot_H_num + H_num;

            int E_num = 0;
			words.clear();
		    NSPutils::split(topnss,words,"CH");
			for (auto &word : words){
				E_num = E_num + word.size();
			}
			tot_E_num = tot_E_num + E_num;

	        int C_num = 0;
			words.clear();
		    NSPutils::split(topnss,words,"EH");
			for (auto &word : words){
				C_num = C_num + word.size();
			}
			tot_C_num = tot_C_num + C_num;
			//std::cout << "num of C" << C_num << std::endl;
		}
	}
	//std::cout << "tot num of H  " << tot_H_num <<"tot num of E  "  << tot_E_num <<"tot num of C  " << tot_C_num << std::endl;
	av_H_num = (double)tot_H_num / (topn+1);
	std::cout << "num of H  " << av_H_num << std::endl;
	if (start_H_num==0) { start_H_num = 1; }
	std::cout << "Percentage change of H  " << (double)(av_H_num-start_H_num)/ start_H_num << std::endl;
	av_E_num = (double)tot_E_num / (topn+1);
	std::cout << "num of E  " << av_E_num << std::endl;
	if (start_E_num==0) { start_E_num = 1; }
	std::cout << "Percentage change of E  " << (double)(av_E_num - start_E_num) / start_E_num << std::endl;
	av_C_num = (double)tot_C_num / (topn + 1);
	std::cout << "num of C  " << av_C_num << std::endl;
	if (start_C_num==0) { start_C_num = 1; }
	std::cout << "Percentage change of C  " << (double)(av_C_num - start_C_num) / start_C_num << std::endl;
}
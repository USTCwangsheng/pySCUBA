/*
 * readrestraints.cpp
 *
 *  Created on: 2021年3月30日
 *      Author: hyiu
 */

#include "sampling/restraints.h"
#include "sampling/rmsdrestraintreader.h"
using namespace NSPsampling;
void readcontactrestraints(NSPintrct::IntrctMol &imol, NSPdataio::InputLines &input,
		int &lidx, std::vector<std::shared_ptr<ResTraint>> &results,int switchkey){
	 NSPintrct::BlkSelector blksa,blksb;
	while (lidx < input.size()) {
			auto & words = input[lidx++];
			if (words[0].substr(0, 3) == "END")
				break;
			if(words[0]=="GroupAResidues"){
				std::string selectedblks;
						for (int i = 1; i < words.size() - 1; ++i) {
							selectedblks = selectedblks + words[i] + " ";
						}
						selectedblks = selectedblks + words.back();
						blksa = NSPintrct::BlkSelector(selectedblks);
						continue;
			} else if(words[0]=="GroupBResidues"){
					std::string selectedblks;
							for (int i = 1; i < words.size() - 1; ++i) {
								selectedblks = selectedblks + words[i] + " ";
							}
							selectedblks = selectedblks + words.back();
							blksb = NSPintrct::BlkSelector(selectedblks);
			} else {
				int nc_0=std::stoi(words[0]);
				double kres=std::stod(words[1]);
				double resgdmin = 5;
				double resgdsmall = 9;
				double resgdoff = 20;
				if (words.size() == 5) {
					resgdmin = std::stod(words[2]);
					resgdsmall = std::stod(words[3]);
					resgdoff = std::stod(words[4]);
				}
				std::shared_ptr<ResTraint>  crest=std::shared_ptr<ResTraint>(
						new ContactRestraint(imol,blksa,blksb,(double) nc_0,kres,resgdmin/10,resgdsmall/10,resgdoff/10));
				crest->switchable() = switchkey;
				crest->restriantname() = "Contact";
				results.push_back(crest);
			}
	}
}
void readscaleinteraction(NSPintrct::IntrctMol &imol, NSPdataio::InputLines &input,
		int &lidx){
	   auto words=input[lidx++];
	   double scale=std::stod(words[0]);
		std::string selectedblks;
		for (int i = 1; i < words.size() - 1; ++i) {
			selectedblks = selectedblks + words[i] + " ";
		}
		selectedblks = selectedblks + words.back();
		auto blks = NSPintrct::BlkSelector(selectedblks).selectedblks;
		for(auto &cm:blks){
			for (auto &r:cm.second){
				auto &blk=imol.getblck(NSPdstl::Idx2D(cm.first,r));
				blk.weight_mc=scale;
			}
		}
}
void readhelixrestraints(NSPintrct::IntrctMol &imol, NSPdataio::InputLines &input,
		int &lidx, std::vector<std::shared_ptr<ResTraint>> &results){
	    auto words=input[lidx++];
	    int nhelices=words.size()/4;
	    int widx=0;
		std::string restriantname{ "Helix" };
	    for (int i=0;i<nhelices;++i){
	    	int chain=std::stoi(words[widx++]);
	    	int start=std::stoi(words[widx++]);
	    	int last=std::stoi(words[widx++]);
	    	double kres=std::stod(words[widx++]);
	    	std::shared_ptr<ResTraint>  hrest=std::shared_ptr<ResTraint>(new HelixRestraint(&imol));
	    	for (int res=start; res<=(last-4); ++res){
	    		auto &blk1= imol.getblck(NSPdstl::Idx2D(chain,res));
	    	    int a1=blk1.aoffset+blk1.natoms()-1;
	   		   int a2=imol.getblck(NSPdstl::Idx2D(chain,res+4)).aoffset;
	   		   ((HelixRestraint*) (hrest.get()))->hbrestraints.push_back(DisRestraint(a1,a2,0.36,0.4,kres));
	    	}
			hrest->switchable() = 0;
			hrest->restriantname() = restriantname;
	    	results.push_back(hrest);
	    }
}
void readrmsdrestraints(NSPintrct::IntrctMol &imol,
		NSPdataio::InputLines &input, int &lidx,
		std::vector<std::shared_ptr<ResTraint>> &results, int& switchkey) {
	RMSDRestraintReader reader;
	reader.readterms(input, lidx);
	reader.addrestraints(imol, results, switchkey);
}
void readrgrestraint(NSPintrct::IntrctMol &imol, NSPdataio::InputLines &input,
		int &lidx, std::vector<std::shared_ptr<ResTraint>> &results,int & switchkey) {
	auto & words = input[lidx++];
	double krg = std::stod(words[0]);
	double rg0 = std::stod(words[1]);
	std::string restriantname{ "RG" };
	std::map<int, std::set<int>> rgblks;
	if (words.size() > 2) {
		if (words[2] == "Rg_residues") {
			std::string selectedblks;
			for (int i = 3; i < words.size() - 1; ++i) {
				selectedblks = selectedblks + words[i] + " ";
			}
			selectedblks = selectedblks + words.back();
			rgblks = NSPintrct::BlkSelector(selectedblks).selectedblks;
		}
		if (krg > 0) {
			std::vector<double> wrg;
			std::vector<int> rgatoms = imol.atomindices(rgblks);
			if (!rgatoms.empty()) {
				wrg.resize(imol.natoms(), 0.0);
				for (auto a : rgatoms)
					wrg[a] = 1.0;
			}
			std::shared_ptr<ResTraint> rgrest = std::shared_ptr<ResTraint>(
					new RgRestraint(rg0 * A2NM, krg * KBT, wrg));
			rgrest->switchable() = switchkey;
			rgrest->restriantname() = restriantname;
			results.push_back(rgrest);
		}
	}
}
std::vector<std::shared_ptr<ResTraint> > NSPsampling::readrestraints(
		NSPintrct::IntrctMol &imol, const std::string & resfile) {
	std::vector<std::shared_ptr<ResTraint>> results;
	NSPdataio::InputLines input;
	input.init(resfile, '#');
	int lidx = 0;
	int switchkey = 0;
	if (input.size() == 0) {
		std::cout << "undefined input in " << resfile << std::endl;
		exit(1);
	}
	while (lidx < input.size()) {
		auto& words = input[lidx++];
		if (words[0] == "RMSDRestraints") {
			switchkey = 0;
//			if (std::stod(words[1]) == 1 || std::stod(words[1]) == 0) { switchkey = std::stod(words[1]); }
			if(words.size()>1){
			if(std::stod(words[1]) == 1 || std::stod(words[1]) == 0){ switchkey = std::stod(words[1]); }
			else { std::cout << "unvalid input in RMSDRestraints: " << std::stod(words[1]) << std::endl;
			exit(1); }
			}
			std::cout << "rmsdswitchkey" << switchkey << std::endl;
			readrmsdrestraints(imol, input, lidx, results, switchkey);
		}
		else if (words[0] == "RgRestraint") {
			switchkey = 0;
//			if(std::stod(words[1])==1|| std::stod(words[1]) == 0){ switchkey = std::stod(words[1]); }	
			if (words.size() > 1) {
				if (std::stod(words[1]) == 1 || std::stod(words[1]) == 0) { switchkey = std::stod(words[1]); }
				else {
					std::cout << "unvalid input in RgRestraint: " << std::stod(words[1]) << std::endl;
					exit(1);
				}
			}
			std::cout << "rgrestraintswitchkey" << switchkey << std::endl;
			readrgrestraint(imol, input, lidx, results, switchkey);
			words = input[lidx++];
			assert(words[0].substr(0, 3) == "END");
		}
		else if (words[0] == "HelixRestraints") {
			switchkey = 0;
			readhelixrestraints(imol, input, lidx, results);
			words = input[lidx++];
			assert(words[0].substr(0, 3) == "END");
		}
		else if (words[0] == "ContactRestraints") {
			switchkey = 0;
//			if (std::stod(words[1]) == 1 || std::stod(words[1]) == 0) { switchkey = std::stod(words[1]); }
			if (words.size() > 1) {
				if (std::stod(words[1]) == 1 || std::stod(words[1]) == 0) { switchkey = std::stod(words[1]); }
				else {
					std::cout << "unvalid input in ContactRestraints: " << std::stod(words[1]) << std::endl;
					exit(1);
				}
			}
			std::cout << "betaswitchkey" << switchkey << std::endl;
			readcontactrestraints(imol, input, lidx, results, switchkey);
		}
		else if (words[0] == "ScaleInteraction") {
			readscaleinteraction(imol, input, lidx);
			words = input[lidx++];
			assert(words[0].substr(0, 3) == "END");
		}
		else {
			std::cout << "undefined input in " << resfile<<". Unsupported key-value found:"<<words[0] << std::endl;
			exit(1);
		}
	}
	return results;
}

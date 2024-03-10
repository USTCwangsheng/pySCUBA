/*
 * loopsamplor.cpp
 *
 *  Created on: 2019年12月24日
 *      Author: hyiu
 */
#include "sampling/loopsampler.h"
#include "dstl/sortindex.h"
#include "iblock/intrctmol.h"
#include "sampling/runsddefaults.h"

using namespace NSPintrct;
using namespace NSPsampling;

BlkSelector::SelectedBlks LoopSampler::inloopblks(
		const std::vector<FlankingSite> &fcites){
	BlkSelector::SelectedBlks sel;
	for(auto &l:fcites){
		sel[l.first.first]=std::set<int>();
	}
	for(auto &l:fcites){
		auto & aset=sel.at(l.first.first);
		for(int i=l.first.second+1; i<l.second.second;++i){
			aset.insert(i);
		}
	}
	return sel;
}
BlkSelector::SelectedBlks LoopSampler::loopblks(const std::vector<std::pair<int, int>>& lranges) {
	BlkSelector::SelectedBlks sel;
	// for(int n=0;n<lranges.size();n++){
	//     sel[n]=std::set<int>();
	// }
	sel[0] = std::set<int>();
	int count = 0;
	for (auto& l : lranges) {
		for (int i = l.first;i < l.second;i++) sel[0].insert(i);
		count++;
	}
	return sel;
}
std::vector<LoopSampler::FlankingSite> NSPsampling::reindexsites(const std::vector<LoopSampler::FlankingSite> & oldsites,
		const std::vector<int> & looplengths){
		std::vector<int> dl;
		for(int i=0;i<oldsites.size();++i){
			dl.push_back(looplengths[i]-(oldsites[i].second.second-oldsites[i].first.second-1));
		}
		std::vector<LoopSampler::FlankingSite> res=oldsites;
		for(int i=0;i<oldsites.size();++i){
			auto & ns=res[i];
			int dposi=0;
			for(int j=0;j<oldsites.size();++j){
				if(oldsites[j].first.first != ns.first.first) continue;
				if(oldsites[j].first.second>=oldsites[i].first.second) continue;
				dposi+=dl[j];
			}
			ns.first.second +=dposi;
			ns.second.second =ns.first.second+looplengths[i]+1;
		}
		return res;
}


std::pair<std::vector<std::pair<MolModeler::SiteIdx, MolModeler::SiteIdx>>, std::vector<std::pair<std::pair<int,int>,int>>> LoopSampler::findloops(NSPintrct::IntrctMol& startmol,const std::vector<FlankingSite>& flankingsites,
	std::vector<int>& looplengths, const LpSamplerPara& lpspara,
	const NSPintrct::IntrctPara& ipara,
	const SDRunPara& sdpara) {
	std::vector<FlankingSite> findflankingsites= flankingsites;
	std::vector<int> findlooplengths;
	std::vector<std::pair<std::pair<int, int>, int>>lrshiftandlengths;
	MolModeler modeler0, modelerfind;//操作0，不操作find
	std::vector<int> start_end_length;
	int looplength = 0;
	std::vector<int> length3; //给每个起始位置固定的loop找3个长度
	int ntest = 0;
	//MolModeler modelerfind;
	modeler0.checkclash(true);
	modelerfind.checkclash(true);
	//NSPdstl::RandomEngine<> rng;
	//rng.init(33); //seed33
	

	int nloops = flankingsites.size();
	std::vector<int> loopindex = NSPdstl::sortindex(flankingsites);//升序排列的loop次序《1，2，3，4》
	imol_ = std::shared_ptr<IntrctMol>(new IntrctMol(startmol));
	modeler0.settarget(*imol_);
	modelerfind.settarget(*imol_);
	std::cout << "Each loop tries 4 starting and ending positions, and each position tries 3 lengths" << std::endl;
	std::cout << "FindLoopTimes = "<<lpspara.findlooptimes << std::endl;
	for (int i = nloops - 1;i >= 0;--i) { //从后往前算loop
		//每段loop设置端点设置4种情况：1原来，2往前加一，3往后加一，4前后各加一
		int idx = loopindex[i];  //loop的次序，从后往前取《4，3，2，1》
		int cachesitefirst = 0;
		int cachesitesecond = 0;
		int templeftshift=0;
		int temprightshift=0;
		int leftshift = 0;
		int rightshift = 0;

		double sigeave_loop = -1000;
		double lastsigetot_loop = 10000;
		double lastsigeave_loop = 10000;
		//给定loop端点和长度的情况
		if (looplengths[i] > 2) {
			//端点不变，长度不变
			looplength=looplengths[i];
			lrshiftandlengths.push_back(std::make_pair(std::make_pair(0,0), looplength));
			std::cout << "use predefine loop:start,end,length" << findflankingsites[idx].first.second << "," << findflankingsites[idx].second.second << "," << looplengths[i] << std::endl;
			continue;
		}
		for (int s = 0;s < 4;s++) {
			if (s == 0) { findflankingsites[idx].first.second = flankingsites[idx].first.second; findflankingsites[idx].second.second = flankingsites[idx].second.second;templeftshift = 0;temprightshift = 0; }
			if (s == 1) { findflankingsites[idx].first.second = flankingsites[idx].first.second - 1; findflankingsites[idx].second.second = flankingsites[idx].second.second;templeftshift = -1;temprightshift = 0; }//左减一
			if (s == 2) { findflankingsites[idx].first.second = flankingsites[idx].first.second; findflankingsites[idx].second.second = flankingsites[idx].second.second + 1;templeftshift = 0;temprightshift = 1; }//右加一
			if (s == 3) {
				findflankingsites[idx].first.second = flankingsites[idx].first.second - 1;
				findflankingsites[idx].second.second = flankingsites[idx].second.second + 1;templeftshift = -1;temprightshift = 1;
			}//左右-+1
			auto confmakerfind = modelerfind.getloopconfmaker(findflankingsites[idx].first, findflankingsites[idx].second);
			//auto confmaker = modeler0.getloopconfmaker(findflankingsites[idx].first, findflankingsites[idx].second);
			auto confmaker = modelerfind.getloopconfmaker(findflankingsites[idx].first, findflankingsites[idx].second);
			std::cout << "test loop sites " << findflankingsites[idx].first.second << ", " << findflankingsites[idx].second.second << std::endl;
			int cid = findflankingsites[idx].first.first;
			int rstart = findflankingsites[idx].first.second - 1;//loop起始位置左移一位-1
			int rend = findflankingsites[idx].second.second + 1;
			int shift = 2;//多算的残基
			if (s == 1) {rstart = rstart + 1;shift = 1;}
			if (s == 2) { rend = rend - 1;shift = 1; }
			if (s == 3) { rstart = rstart + 1;rend = rend - 1;shift = 0; }
			std::vector<IntrctBlck> newloop;
			length3.clear();
			for (int i = 3;i < 20;i++) {
				//newloop = modeler0.newloop(confmaker, std::vector<std::string>(i, "GLY"), 0, 100);
				newloop = modelerfind.newloop(confmakerfind, std::vector<std::string>(i, "GLY"), 0, 100);
				std::cout << "test loop length " << i << std::endl;
				if (newloop.empty()) { continue; }
				else {
					length3.push_back(i);
					if (length3.size() == 3) break;
				}
			}
			if (!length3.empty()) {
				for (auto& nlen : length3) {
					int times = 0;
					double sigetot_loop = 0;
					BlkSelector::SelectedBlks mclr;
					std::vector<std::pair<int, int>> lr;
					int loopend = rstart + nlen + 3 - 1;
					if (s == 1) { loopend = loopend - 1; };
					if (s == 2) { loopend = loopend - 1; };
					if (s == 3) { loopend = loopend - 2; };
					lr.push_back(std::make_pair(rstart + 1, loopend));//确定可以活动的残基，loop起点
					mclr = loopblks(lr);
					while (times < lpspara.findlooptimes) {
						std::cout << " in loop: " << i << " loop length: " << nlen << " times: " << times << std::endl;
						modeler0 = modelerfind;
						newloop = modeler0.newloop(confmaker, std::vector<std::string>(nlen, "GLY"), 0, 200);
						if (!newloop.empty()) {
							times++;
							auto imol = std::shared_ptr<IntrctMol>(new IntrctMol(
								modeler0.mkimol_newloop(findflankingsites[idx].first, findflankingsites[idx].second, newloop)));
							//std::cout << "mkimol_newloop "  << std::endl;
							modeler0.settarget(*imol);
							//int loopend = rstart + nlen + 3 - 1;
							//if (s == 1) { loopend = loopend - 1; };
							//if (s == 2) { loopend = loopend - 1; };
							//if (s == 3) { loopend = loopend - 2; };
							//lr.push_back(std::make_pair(rstart + 1, loopend));//确定可以活动的残基，loop起点
							//mclr = loopblks(lr);
							imol->specifyactiveblks(mclr, typename NSPintrct::BlkSelector::SelectedBlks());
							imol->my_intrctpara = std::shared_ptr<IntrctPara>(new IntrctPara(mkrescaledparam(ipara)));
							optimizer_.init(sdpara, imol);
							const std::vector<double>& initcrd = imol->getcrds_all_d();
							optimizer_.reinitstate(initcrd, 0);
							SDCallBackOptSD callbck(mypara_);
							optimizer_.runsteps(mypara_.maxsdsteps, callbck);
							for (int r = rstart + 1; r <= loopend; ++r) { //int r = rstart; r <= rstart+nlen+3; ++r 现在往左右加一
								//std::cout << "residues: " << r << std::endl;
								for (int t = 0;t < NSPintrct::IntrctBlck::ENESIZE;++t) {
									//std::cout << "energies[t] " << imol->getblck(NSPdstl::Idx2D(cid, r)).energies[t] << std::endl;
									sigetot_loop += imol->getblck(NSPdstl::Idx2D(cid, r)).energies[t]; //ene of all residues in loops-2+2 *100
									//std::cout << " in loop  "<<i<<"sigetot_loop" << sigetot_loop << std::endl;
								}
							}

						}
						else { std::cout << "find length but cant rebuild in findloop" << std::endl;	continue; }
					}

					sigeave_loop = sigetot_loop / (lpspara.findlooptimes * (nlen+ shift));
					//std::cout << "count nlen+shift: " << nlen+ shift << std::endl;
					std::cout << "leng :" << nlen << "sigeave_loop :" << sigeave_loop << std::endl;
					if (sigeave_loop < lastsigeave_loop) {
						std::cout << "lower ene  sigeave_loop: " << sigeave_loop << " < lastsigeave_loop : " << lastsigeave_loop << " nlen = " << nlen << std::endl;
						lastsigeave_loop = sigeave_loop;
						looplength = nlen;
						cachesitefirst = findflankingsites[idx].first.second;
						cachesitesecond = findflankingsites[idx].second.second;
						leftshift = templeftshift;
						rightshift = temprightshift;
						std::cout << "cachesite " << cachesitefirst << " , " << cachesitesecond << std::endl;
					}//在（小于等于）三个不同长度内取最低能量
				}
			}
		}
		findflankingsites[idx].first.second = cachesitefirst;
		findflankingsites[idx].second.second = cachesitesecond;
		//findlooplengths.push_back(looplength);
		auto lrshift=std::make_pair(leftshift, rightshift);
		lrshiftandlengths.push_back(std::make_pair(lrshift, looplength));
		std::cout << "find best loop length: " << looplength << " in sites" << findflankingsites[idx].first.second << ", " << findflankingsites[idx].second.second << std::endl;
	}
	reverse(lrshiftandlengths.begin(), lrshiftandlengths.end());
	auto sitesandlength = std::make_pair(findflankingsites, lrshiftandlengths);
	return sitesandlength;
}

std::pair<std::vector<std::pair<MolModeler::SiteIdx, MolModeler::SiteIdx>>, std::vector<std::pair<std::pair<int, int>, int>>> LoopSampler::findshortloops(NSPintrct::IntrctMol& startmol, const std::vector<FlankingSite>& flankingsites,
	std::vector<int>& looplengths, const LpSamplerPara& lpspara,
	const NSPintrct::IntrctPara& ipara,
	const SDRunPara& sdpara) {
	std::vector<FlankingSite> findflankingsites = flankingsites;//此容器计算从0开始，我们的loop从1开始
	std::vector<int> findlooplengths;
	std::vector<std::pair<std::pair<int, int>, int>>lrshiftandlengths;
	MolModeler modeler0, modelerfind;//操作0，不操作find
	std::vector<int> start_end_length;
	int looplength = 0;
	std::vector<int> length3; //给每个起始位置固定的loop找3个长度
	int shortlength;
	int ntest = 0;
	//MolModeler modelerfind;
	modeler0.checkclash(true);
	modelerfind.checkclash(true);
	//NSPdstl::RandomEngine<> rng;
	//rng.init(33); //seed33


	int nloops = flankingsites.size();
	std::vector<int> loopindex = NSPdstl::sortindex(flankingsites);//升序排列的loop次序《1，2，3，4》
	imol_ = std::shared_ptr<IntrctMol>(new IntrctMol(startmol));
	modeler0.settarget(*imol_);
	modelerfind.settarget(*imol_);
	std::cout << "Try 4 starting and ending positions for each loop, try the shortest length for each position" << std::endl;
	std::cout << "FindLoopTimes = " << lpspara.findlooptimes << std::endl;
	for (int i = nloops - 1;i >= 0;--i) { //从后往前算loop
		//每段loop设置端点设置4种情况：1原来，2往前加一，3往后加一，4前后各加一
		int idx = loopindex[i];  //loop的次序，从后往前取《4，3，2，1》
		int cachesitefirst = 0;
		int cachesitesecond = 0;
		int templeftshift = 0;
		int temprightshift = 0;
		int leftshift = 0;
		int rightshift = 0;

		double sigeave_loop = -1000;
		double lastsigetot_loop = 10000;
		double lastsigeave_loop = 10000;
		//给定loop端点和长度的情况
		if (looplengths[i] > 2) {
			//端点不变，长度不变
			looplength = looplengths[i];
			lrshiftandlengths.push_back(std::make_pair(std::make_pair(0, 0), looplength));
			std::cout << "use predefine loop:start,end,length" << findflankingsites[idx].first.second << "," << findflankingsites[idx].second.second << "," << looplengths[i] << std::endl;
			continue;
		}
		for (int s = 0;s < 4;s++) {
			if (s == 0) { findflankingsites[idx].first.second = flankingsites[idx].first.second; findflankingsites[idx].second.second = flankingsites[idx].second.second;templeftshift = 0;temprightshift = 0;}
			if (s == 1) { findflankingsites[idx].first.second = flankingsites[idx].first.second - 1; findflankingsites[idx].second.second = flankingsites[idx].second.second;templeftshift = -1;temprightshift = 0;}//左减一
			if (s == 2) { findflankingsites[idx].first.second = flankingsites[idx].first.second; findflankingsites[idx].second.second = flankingsites[idx].second.second + 1;templeftshift = 0;temprightshift = 1; }//右加一
			if (s == 3) {
				findflankingsites[idx].first.second = flankingsites[idx].first.second - 1;
				findflankingsites[idx].second.second = flankingsites[idx].second.second + 1;templeftshift = -1;temprightshift = 1;
			}//左右-+1
			auto confmakerfind = modelerfind.getloopconfmaker(findflankingsites[idx].first, findflankingsites[idx].second);
			//auto confmaker = modeler0.getloopconfmaker(findflankingsites[idx].first, findflankingsites[idx].second);
			auto confmaker = modelerfind.getloopconfmaker(findflankingsites[idx].first, findflankingsites[idx].second);
			std::cout << "test loop sites " << findflankingsites[idx].first.second << ", " << findflankingsites[idx].second.second << std::endl;
			int cid = findflankingsites[idx].first.first;
			int rstart = findflankingsites[idx].first.second - 1;//loop起始位置左移一位-1，用于计算平均能量
			int rend = findflankingsites[idx].second.second + 1;//loop终止位置右移一位+1
			int shift = 2;//多算的残基
			if (s == 1) { rstart = rstart + 1;shift = 1; }
			if (s == 2) { rend = rend - 1;shift = 1; }
			if (s == 3) { rstart = rstart + 1;rend = rend - 1;shift = 0; }
			std::vector<IntrctBlck> newloop;
			shortlength=0;
			for (int i = 3;i < 20;i++) {
				//newloop = modeler0.newloop(confmaker, std::vector<std::string>(i, "GLY"), 0, 100);
				newloop = modelerfind.newloop(confmakerfind, std::vector<std::string>(i, "GLY"), 0, 100);
				std::cout << "test loop length " << i << std::endl;
				if (newloop.empty()) { continue; }
				else {
					shortlength = i;
					break;
				}
			}
			if (shortlength!=0) {
					int times = 0;
					double sigetot_loop = 0;
					BlkSelector::SelectedBlks mclr;
					std::vector<std::pair<int, int>> lr;
					int loopend = rstart + shortlength + 3 - 1;
					if (s == 1) { loopend = loopend - 1; };
					if (s == 2) { loopend = loopend - 1; };
					if (s == 3) { loopend = loopend - 2; };
					lr.push_back(std::make_pair(rstart + 1, loopend));//确定可以活动的残基，loop起点,这里每次只有一段loop
					mclr = loopblks(lr);
					while (times < lpspara.findlooptimes) {
						std::cout << " in loop: " << i << " loop length: " << shortlength << " times: " << times << std::endl;
						modeler0 = modelerfind;
						newloop = modeler0.newloop(confmaker, std::vector<std::string>(shortlength, "GLY"), 0, 200);
						if (!newloop.empty()) {
							times++;
							auto imol = std::shared_ptr<IntrctMol>(new IntrctMol(
								modeler0.mkimol_newloop(findflankingsites[idx].first, findflankingsites[idx].second, newloop)));
							//std::cout << "mkimol_newloop "  << std::endl;
							modeler0.settarget(*imol);
							//int loopend = rstart + shortlength + 3 - 1;
							//if (s == 1) { loopend = loopend - 1; };
							//if (s == 2) { loopend = loopend - 1; };
							//if (s == 3) { loopend = loopend - 2; };
							//lr.push_back(std::make_pair(rstart + 1, loopend));//确定可以活动的残基，loop起点
							//mclr = loopblks(lr);
							imol->specifyactiveblks(mclr, typename NSPintrct::BlkSelector::SelectedBlks());
							imol->my_intrctpara = std::shared_ptr<IntrctPara>(new IntrctPara(mkrescaledparam(ipara)));
							optimizer_.init(sdpara, imol);
							const std::vector<double>& initcrd = imol->getcrds_all_d();
							optimizer_.reinitstate(initcrd, 0);
							SDCallBackOptSD callbck(mypara_);
							optimizer_.runsteps(mypara_.maxsdsteps, callbck);
							for (int r = rstart + 1; r <= loopend; ++r) { //int r = rstart; r <= rstart+nlen+3; ++r 现在往左右加一
								//std::cout << "residues: " << r << std::endl;
								for (int t = 0;t < NSPintrct::IntrctBlck::ENESIZE;++t) {
									//std::cout << "energies[t] " << imol->getblck(NSPdstl::Idx2D(cid, r)).energies[t] << std::endl;
									sigetot_loop += imol->getblck(NSPdstl::Idx2D(cid, r)).energies[t]; //ene of all residues in loops-2+2 *100
									//std::cout << " in loop  "<<i<<"sigetot_loop" << sigetot_loop << std::endl;
								}
							}

						}
						else { std::cout << "find length but cant rebuild in findloop" << std::endl;	continue; }
					}

					sigeave_loop = sigetot_loop / (lpspara.findlooptimes * (shortlength + shift));
					//std::cout << "count shortlength+shift: " << shortlength+ shift << std::endl;
					std::cout << "leng :" << shortlength << "sigeave_loop :" << sigeave_loop << std::endl;
					if (sigeave_loop < lastsigeave_loop) {
						std::cout << "lower ene  sigeave_loop: " << sigeave_loop << " < lastsigeave_loop : " << lastsigeave_loop << " nlen = " << shortlength << std::endl;
						lastsigeave_loop = sigeave_loop;
						looplength = shortlength;
						cachesitefirst = findflankingsites[idx].first.second;
						cachesitesecond = findflankingsites[idx].second.second;
						leftshift = templeftshift;
						rightshift = temprightshift;
						std::cout << "cachesite " << cachesitefirst << " , " << cachesitesecond << std::endl;
					}//在（小于等于）三个不同长度内取最低能量
			}
		}
		findflankingsites[idx].first.second = cachesitefirst;
		findflankingsites[idx].second.second = cachesitesecond;
		//findlooplengths.push_back(looplength);
		auto lrshift = std::make_pair(leftshift, rightshift);
		lrshiftandlengths.push_back(std::make_pair(lrshift, looplength));
		std::cout << "find best loop length: " << looplength << " in sites" << findflankingsites[idx].first.second << ", " << findflankingsites[idx].second.second << std::endl;
	}
	reverse(lrshiftandlengths.begin(), lrshiftandlengths.end());
	auto sitesandlength = std::make_pair(findflankingsites, lrshiftandlengths);
	return sitesandlength;
}

std::pair<std::vector<std::pair<MolModeler::SiteIdx, MolModeler::SiteIdx>>, std::vector<std::pair<std::pair<int, int>, int>>> LoopSampler::fixloops(NSPintrct::IntrctMol& startmol, const std::vector<FlankingSite>& flankingsites,
	std::vector<int>& looplengths, const LpSamplerPara& lpspara,
	const NSPintrct::IntrctPara& ipara,
	const SDRunPara& sdpara) {
	std::vector<FlankingSite> findflankingsites = flankingsites;
	std::vector<int> findlooplengths;
	std::vector<std::pair<std::pair<int, int>, int>>lrshiftandlengths;
	MolModeler modeler0, modelerfind;//操作0，不操作find
	std::vector<int> start_end_length;
	int looplength = 0;
	std::vector<int> length3; //给每个起始位置固定的loop找3个长度
	int ntest = 0;
	modeler0.checkclash(true);
	modelerfind.checkclash(true);

	int nloops = flankingsites.size();
	std::vector<int> loopindex = NSPdstl::sortindex(flankingsites);//升序排列的loop次序《1，2，3，4》
	imol_ = std::shared_ptr<IntrctMol>(new IntrctMol(startmol));
	modeler0.settarget(*imol_);
	modelerfind.settarget(*imol_);
	std::cout << "Try fixed starting and ending positions for each loop, and each position tries 3 lengths" << std::endl;
	std::cout << "FindLoopTimes = " << lpspara.findlooptimes << std::endl;
	for (int i = nloops - 1; i >= 0; --i) { //从后往前算loop
		//每段loop设置端点设置4种情况：1原来，2往前加一，3往后加一，4前后各加一
		int idx = loopindex[i];  //loop的次序，从后往前取《4，3，2，1》
		int leftshift = 0;
		int rightshift = 0;
		double sigeave_loop = -1000;
		double lastsigetot_loop = 10000;
		double lastsigeave_loop = 10000;
		//给定loop端点和长度的情况
		if (looplengths[i] > 2) {
			//端点不变，长度不变
			looplength = looplengths[i];
			lrshiftandlengths.push_back(std::make_pair(std::make_pair(0, 0), looplength));
			std::cout << "use predefine loop:start,end,length" << findflankingsites[idx].first.second << "," << findflankingsites[idx].second.second << "," << looplengths[i] << std::endl;
			continue;
		}
			auto confmakerfind = modelerfind.getloopconfmaker(findflankingsites[idx].first, findflankingsites[idx].second);
			//auto confmaker = modeler0.getloopconfmaker(findflankingsites[idx].first, findflankingsites[idx].second);
			auto confmaker = modelerfind.getloopconfmaker(findflankingsites[idx].first, findflankingsites[idx].second);
			std::cout << "test loop sites " << findflankingsites[idx].first.second << ", " << findflankingsites[idx].second.second << std::endl;
			int cid = findflankingsites[idx].first.first;
			int rstart = findflankingsites[idx].first.second - 1;//loop起始位置往前多算一位-1，用于计算平均能量
			int rend = findflankingsites[idx].second.second + 1;//往后多算一位
			int shift = 2;//多算的残基

			std::vector<IntrctBlck> newloop;
			length3.clear();
			for (int i = 3; i < 20; i++) {
				//找loop长度
				newloop = modelerfind.newloop(confmakerfind, std::vector<std::string>(i, "GLY"), 0, 100);
				std::cout << "test loop length " << i << std::endl;
				if (newloop.empty()) { continue; }
				else {
					length3.push_back(i);
					if (length3.size() == 3) break;
				}
			}
			if (!length3.empty()) {
				for (auto& nlen : length3) {
					int times = 0;
					double sigetot_loop = 0;
					BlkSelector::SelectedBlks mclr;
					std::vector<std::pair<int, int>> lr;
					int loopend = rstart + nlen + 3 - 1;
					lr.push_back(std::make_pair(rstart + 1, loopend));//确定可以活动的残基，loop起点
					mclr = loopblks(lr);
					while (times < lpspara.findlooptimes) {
						std::cout << " in loop: " << i << " loop length: " << nlen << " times: " << times << std::endl;
						modeler0 = modelerfind;
						newloop = modeler0.newloop(confmaker, std::vector<std::string>(nlen, "GLY"), 0, 200);
						if (!newloop.empty()) {
							times++;
							auto imol = std::shared_ptr<IntrctMol>(new IntrctMol(
								modeler0.mkimol_newloop(findflankingsites[idx].first, findflankingsites[idx].second, newloop)));
							//std::cout << "mkimol_newloop "  << std::endl;
							modeler0.settarget(*imol);
							imol->specifyactiveblks(mclr, typename NSPintrct::BlkSelector::SelectedBlks());
							imol->my_intrctpara = std::shared_ptr<IntrctPara>(new IntrctPara(mkrescaledparam(ipara)));
							optimizer_.init(sdpara, imol);
							const std::vector<double>& initcrd = imol->getcrds_all_d();
							optimizer_.reinitstate(initcrd, 0);
							SDCallBackOptSD callbck(mypara_);
							optimizer_.runsteps(mypara_.maxsdsteps, callbck);
							for (int r = rstart + 1; r <= loopend; ++r) { //int r = rstart; r <= rstart+nlen+3; ++r 现在往左右加一
								//std::cout << "residues: " << r << std::endl;
								for (int t = 0; t < NSPintrct::IntrctBlck::ENESIZE; ++t) {
									//std::cout << "energies[t] " << imol->getblck(NSPdstl::Idx2D(cid, r)).energies[t] << std::endl;
									sigetot_loop += imol->getblck(NSPdstl::Idx2D(cid, r)).energies[t]; //ene of all residues in loops-2+2 *100
									//std::cout << " in loop  "<<i<<"sigetot_loop" << sigetot_loop << std::endl;
								}
							}

						}
						else { std::cout << "find length but cant rebuild in findloop" << std::endl;	continue; }
					}

					sigeave_loop = sigetot_loop / (lpspara.findlooptimes * (nlen + shift));
					//std::cout << "count nlen+shift: " << nlen+ shift << std::endl;
					std::cout << "leng :" << nlen << "sigeave_loop :" << sigeave_loop << std::endl;
					if (sigeave_loop < lastsigeave_loop) {
						std::cout << "lower ene  sigeave_loop: " << sigeave_loop << " < lastsigeave_loop : " << lastsigeave_loop << " nlen = " << nlen << std::endl;
						lastsigeave_loop = sigeave_loop;
						looplength = nlen;
					}//在（小于等于）三个不同长度内取最低能量
				}
			}
			else {
				std::cout << "cant find fixshortloop in site: " << findflankingsites[idx].first.second << ", " << findflankingsites[idx].second.second << std::endl;
				exit(1);
			}		
		auto lrshift = std::make_pair(leftshift, rightshift);
		lrshiftandlengths.push_back(std::make_pair(lrshift, looplength));
		std::cout << "find best loop length: " << looplength << " in sites" << findflankingsites[idx].first.second << ", " << findflankingsites[idx].second.second << std::endl;
	}
	reverse(lrshiftandlengths.begin(), lrshiftandlengths.end());
	auto sitesandlength = std::make_pair(findflankingsites, lrshiftandlengths);
	return sitesandlength;
}
std::pair<std::vector<std::pair<MolModeler::SiteIdx, MolModeler::SiteIdx>>, std::vector<std::pair<std::pair<int, int>, int>>> LoopSampler::fixshortloops(NSPintrct::IntrctMol& startmol, const std::vector<FlankingSite>& flankingsites,
	std::vector<int>& looplengths, const LpSamplerPara& lpspara,
	const NSPintrct::IntrctPara& ipara,
	const SDRunPara& sdpara) {
	std::vector<FlankingSite> findflankingsites = flankingsites;
	std::vector<int> findlooplengths;
	std::vector<std::pair<std::pair<int, int>, int>>lrshiftandlengths;
	MolModeler modelerfind;//操作0，不操作find
	std::vector<int> start_end_length;
	int looplength = 0;
	int shortlength;
	modelerfind.checkclash(true);

	int nloops = flankingsites.size();
	std::vector<int> loopindex = NSPdstl::sortindex(flankingsites);//升序排列的loop次序《1，2，3，4》
	imol_ = std::shared_ptr<IntrctMol>(new IntrctMol(startmol));
	modelerfind.settarget(*imol_);
	std::cout << "Try fixed starting and ending positions for each loop, try the shortest length for each position" << std::endl;
	std::cout << "FindLoopTimes =  only 1 time"  << std::endl;
	for (int i = nloops - 1;i >= 0;--i) { //从后往前算loop
		//每段loop设置端点固定
		int idx = loopindex[i];  //loop的次序，从后往前取《4，3，2，1》
		int leftshift = 0;
		int rightshift = 0;

		//给定loop端点和长度的情况
		if (looplengths[i] > 2) {
			//端点不变，长度不变
			looplength = looplengths[i];
			lrshiftandlengths.push_back(std::make_pair(std::make_pair(0, 0), looplength));
			std::cout << "use predefine loop:start,end,length" << findflankingsites[idx].first.second << "," << findflankingsites[idx].second.second << "," << looplengths[i] << std::endl;
			continue;
		}

			auto confmakerfind = modelerfind.getloopconfmaker(findflankingsites[idx].first, findflankingsites[idx].second);
			std::cout << "test loop sites " << findflankingsites[idx].first.second << ", " << findflankingsites[idx].second.second << std::endl;//位置不变
			std::vector<IntrctBlck> newloop;
			shortlength = 0;
			for (int i = 3;i < 20;i++) {
				//找最短loop长度
				newloop = modelerfind.newloop(confmakerfind, std::vector<std::string>(i, "GLY"), 0, 100);
				std::cout << "test loop length " << i << std::endl;
				if (newloop.empty()) { continue; }
				else {
					shortlength = i;
					break;
				}
			}
			if (shortlength != 0) {
				auto lrshift = std::make_pair(leftshift, rightshift);
				lrshiftandlengths.push_back(std::make_pair(lrshift, shortlength));
				std::cout << "find best loop length: " << shortlength << " in sites" << findflankingsites[idx].first.second << ", " << findflankingsites[idx].second.second << std::endl;
			}
			else{ std::cout <<"cant find fixshortloop in site: "<< findflankingsites[idx].first.second << ", " << findflankingsites[idx].second.second << std::endl;
			exit(1);
			}
	}
	reverse(lrshiftandlengths.begin(), lrshiftandlengths.end());
	auto sitesandlength = std::make_pair(findflankingsites, lrshiftandlengths);
	return sitesandlength;
}

bool LoopSampler::setupmol(NSPintrct::IntrctMol& startmol,
	const std::vector<FlankingSite>& flankingsites,
	const std::vector<int>& looplengths) {
	MolModeler modeler0;
	modeler0.checkclash(true);
	int nloops = flankingsites.size();
	std::vector<int> loopindex = NSPdstl::sortindex(flankingsites);
	imol_ = std::shared_ptr<IntrctMol>(new IntrctMol(startmol));
	modeler0.settarget(*imol_);
	for (int i = nloops - 1;i >= 0;--i) {
		int idx = loopindex[i];
		if (flankingsites[idx].second.second - flankingsites[idx].first.second - 1 == looplengths[idx]) {
			for (int r = flankingsites[idx].first.second + 1; r < flankingsites[idx].second.second;++r) {
				modeler0.changeresidue(NSPdstl::Idx2D(flankingsites[idx].first.first, r), "GLY");
			}
			continue;  //do not rebuild loop at beginning if loop length unchanged
		}
		auto confmaker = modeler0.getloopconfmaker(flankingsites[idx].first,
			flankingsites[idx].second);
		auto newloop = modeler0.newloop(confmaker, std::vector<std::string>(looplengths[idx], "GLY"), 99, 1000);
		if (newloop.empty()) {
			std::cout << "Generating loop of length " << looplengths[idx] << "between positions"
				<< flankingsites[idx].first.second << " and " << flankingsites[idx].second.second
				<< " in chain " << flankingsites[idx].first.first << " failed." << std::endl;
			return false;
		}
		imol_ = std::shared_ptr<IntrctMol>(new IntrctMol(
			modeler0.mkimol_newloop(flankingsites[idx].first, flankingsites[idx].second, newloop)));
		modeler0.settarget(*imol_);
	}

	std::cout << "New number of atoms: " << imol_->natoms() << std::endl;
	modeler_.settarget(*imol_);
	reindexedflankingsites_ = reindexsites(flankingsites, looplengths);
	looplengths_ = looplengths;
	loopatoms_.clear();
	for (auto& fs : reindexedflankingsites_) {
		int cid = fs.first.first;
		for (int r = fs.first.second + 1;r < fs.second.second;++r) {
			auto& loopblk = (*imol_)({ cid,r });
			for (int i = 0;i < loopblk.natoms();i++)loopatoms_.push_back(i + loopblk.aoffset);
		}
	}
	for (int i = 0;i < nloops;++i) {
		confmakers_.push_back(modeler_.getloopconfmaker(reindexedflankingsites_[i].first,
			reindexedflankingsites_[i].second));
	}
	return true;
}

	


bool LoopSampler::setupmolfind(NSPintrct::IntrctMol &startmol,
		const std::vector<FlankingSite> & flankingsites,
		std::vector<int> & looplengths){
	MolModeler modeler0;
	modeler0.checkclash(true);
	int nloops=flankingsites.size();
	std::vector<int> loopindex=NSPdstl::sortindex(flankingsites);//升序排列的loop次序《1，2，3，4》
	imol_=std::shared_ptr<IntrctMol>(new IntrctMol(startmol));
	modeler0.settarget(*imol_);
	for (int i=nloops-1;i>=0;--i) { //从后往前算loop
		int idx=loopindex[i];  //loop的次序，从后往前取《4，3，2，1》
		auto confmaker=modeler0.getloopconfmaker(flankingsites[idx].first,flankingsites[idx].second);
		std::vector<IntrctBlck> newloop;
		newloop = modeler0.newloop(confmaker, std::vector<std::string>(looplengths[idx], "GLY"), 0, 1000);
		if (newloop.empty()){
				std::cout << "loop find " << looplengths[idx] << " between positions"
					<< flankingsites[idx].first.second + 1 << " and " << flankingsites[idx].second.second + 1
					<< " in chain " << flankingsites[idx].first.first <<" bad."<< std::endl;
				return false;
		}
		imol_=std::shared_ptr<IntrctMol>(new IntrctMol(
				modeler0.mkimol_newloop(flankingsites[idx].first,flankingsites[idx].second,newloop)));
		modeler0.settarget(*imol_);
	}

	std::cout <<"New number of atoms: " << imol_->natoms() <<std::endl;
	modeler_.settarget(*imol_);
	reindexedflankingsites_=reindexsites(flankingsites,looplengths);
	looplengths_=looplengths;
	loopatoms_.clear();
	for(auto &fs:reindexedflankingsites_){
		int cid=fs.first.first;
		for(int r=fs.first.second+1;r<fs.second.second;++r){
			auto & loopblk=(*imol_)({cid,r});
			for(int i=0;i<loopblk.natoms();i++)loopatoms_.push_back(i+loopblk.aoffset);
		}
	}
	for(int i=0;i<nloops;++i){
		confmakers_.push_back(modeler_.getloopconfmaker(reindexedflankingsites_[i].first,
				reindexedflankingsites_[i].second));
	}
	return true;
}
bool LoopSampler::rebuildloops(const std::vector<int> &lidx){
	for(auto i:lidx){
		std::cout << "rebuilloop: " << i << std::endl;
		auto newloop=modeler_.newloop(confmakers_[i],
				std::vector<std::string>(looplengths_[i],"GLY"),0,200,helixinloop_);
		if (newloop.empty()) {
			std::cout << "cannotrebuildloop " << i << std::endl;
			std::vector<int> r;
			r.push_back(i);
			reuseloops(r);
			std::cout << "reuseloop in rebuild step " << i << std::endl;
		}
		else {
			//test
			//std::cout << "newloop not empty " << i << std::endl;
			modeler_.replaceloop(reindexedflankingsites_[i].first,
				reindexedflankingsites_[i].second, newloop);
			std::cout << "rebuildone: " << i << std::endl;
		}
	}
	std::cout << "readrebuild all loops_done " << std::endl;
	return true;
}
bool LoopSampler::reuseloops(const std::vector<int> &lidy){
	for(auto i:lidy){
		int j = rand() % looppools_[i].loopconfs.size();
		// test_cyx
		std::cout << "LoopPool " << i << " useloop " << j << std::endl;
		for (auto it = looppools_[i].loopconfs.begin(); it != looppools_[i].loopconfs.end(); it++)
		{
			if (j == 0)
			{
				auto newloop=it->second;
				if(newloop.empty()) return false;
				modeler_.replaceloop(reindexedflankingsites_[i].first,
						reindexedflankingsites_[i].second, newloop);
				//test
				//std::cout << "single loop reuselooppools_done" << i << std::endl;
				break;
			}
			j--;
		}
	}
	if(!lidy.empty()) std::cout << "readreuselooppools_done " << std::endl;
	return true;
}
void LoopSampler::optimizeloops(){
	const std::vector<double> &initcrd=imol_->getcrds_all_d();
	optimizer_.reinitstate(initcrd,0);
	SDCallBackOptSD callbck(mypara_);
	optimizer_.runsteps(mypara_.maxsdsteps,callbck);
//	std::vector<NSPgeometry::XYZ>dedx;
// last call for forces
//	imol_->forces_all(NSPgeometry::doublevtoXYZv(optimizer_.state().crd),dedx);
}

std::vector<double> LoopSampler::loopenes(std::vector<std::array<double,NSPintrct::IntrctBlck::ENESIZE>> & energies)
 const {
	int nloops=looplengths_.size();
	std::vector<double> enes(nloops,0.0);
	std::array<double,NSPintrct::IntrctBlck::ENESIZE> etmp;
	for(auto &t:etmp) t=0.0;
	energies=std::vector<std::array<double,NSPintrct::IntrctBlck::ENESIZE>>(nloops,etmp);
	for(int l=0;l<nloops;++l){
		int cid=this->reindexedflankingsites_[l].first.first;
		int rstart=this->reindexedflankingsites_[l].first.second-1;
		if(rstart<0) rstart=0;
		int rend=this->reindexedflankingsites_[l].second.second+1;
		if(rend>=imol_->nresidues(cid)) rend=imol_->nresidues(cid)-1;
		for(int r=rstart; r<=rend; ++r){
			for(int t=0;t<NSPintrct::IntrctBlck::ENESIZE;++t){
				energies[l][t] +=imol_->getblck(NSPdstl::Idx2D(cid,r)).energies[t];
			}
		}
		for(auto t:energies[l]) enes[l]+=t;
	}
	return enes;
}
#include "dataio/parameters.h"
LpSamplerPara::LpSamplerPara(const std::vector<std::string> &controllines){
	   std::map<std::string, double> doublepars {{"ENEAverageDecay",0.99},
	   {"ENEVarStop",10.0},{"MCkBT",20.0}};
		std::map<std::string, std::vector<std::string>> stringvecpars {{"Loops",{}}};
		std::map<std::string, std::vector<double>> doublevecpars {};
		std::map<std::string, std::string> stringpars {{"JobName",""},{"PDBStart",""},{"LoopSearchMode",""},
			{"OutLoopFile",""},{"PMax",""},{"UnchangedCut",""},{"PoolSizeCut",""},{"ReadRMSD",""},{"PhiPsiRMSD",""}};
		std::map<std::string, std::vector<int>> intvecpars {};
		std::map<std::string, int> intpars{{"MaxSDSteps",0},{"MCRandomSeed",1357},
			{"MaxMCCycles",0},{"MCCyclesCut",20},{"PrintParameters",0},{"Verbose",1},
			{"ReconstructNum",0},{"ReadLoopNum",0},{"HelixInLoop",0},{"FindLoopTimes",50}};
	NSPdataio::ParameterSet pset;
	pset.initdefaultkeys(doublepars,stringpars,intpars,doublevecpars,stringvecpars,intvecpars);
	pset.adjustvalues(controllines);
	pset.getval("ENEAverageDecay",&(this->enedecay));
	pset.getval("ENEVarStop",&(this->sdenevarcut));
	pset.getval("MaxSDSteps",&(this->maxsdsteps));
	pset.getval("MaxMCCycles",&(this->maxcycles));
	pset.getval("MCCyclesCut",&(this->MCcyclescut));
	pset.getval("JobName",&(this->jobname));
	pset.getval("PDBStart",&(this->pdbstart));
	pset.getval("LoopSearchMode", &(this->loopsearchmode));
	pset.getval("FindLoopTimes", &(this->findlooptimes));
	pset.getval("OutLoopFile",&(this->outloopfile));
	pset.getval("MCRandomSeed", &(this->mcrandomseed));
	pset.getval("MCkBT", &(this->mckbt));
	pset.getval("PMax", &(this->pmax));
	pset.getval("UnchangedCut", &(this->unchangedcut));
	pset.getval("PoolSizeCut", &(this->poolsizecut));
	pset.getval("ReconstructNum", &(this->reconstructnum));
	pset.getval("ReadLoopNum", &(this->readloopnum));
	pset.getval("ReadRMSD", &(this->readrmsd));
	pset.getval("PhiPsiRMSD", &(this->phipsirmsd));
	std::vector<std::string> loops;
	pset.getval("Loops",&loops);
	assert(loops.size()%4==0);
	int idx=0;
	this->nloops=loops.size()/4;
	for(int i=0;i<this->nloops;i++){
		char cid=loops[idx++][0];
		if(cid=='0') cid=' ';
		std::string left=loops[idx++];
		std::string right=loops[idx++];
		int length=std::stoi(loops[idx++]);
		this->flankingpositions.push_back({{cid,left},{cid,right}});
		std::cout <<"Loop " <<i << " chain " <<cid <<" " <<left <<"- "<<right<< " length " <<length<<std::endl;
//		flankingsites.push_back({{cid,left},{cid,right}});
		looplengths.push_back(length);
	}
	int printpara;
	pset.getval("PrintParameters",&printpara);
	if(printpara !=0){
		auto & oss=NSPintrct::Results::ofstreams(this->jobname,"_par",".dat");
		oss<<"START LoopSamplerPar"<<std::endl;
		pset.readline("PrintParameters = 0");
		pset.printparameters(oss);
		oss<<"END LoopSamplerPar"<<std::endl;
	}
	int verbosein;
	pset.getval("Verbose",&verbosein);
	verbose=true;
	if(verbosein==0) verbose=false;
	int hil;
	pset.getval("HelixInLoop", &hil);
	helixinloop = false;
	if (hil == 1) helixinloop = true;
}

void LoopSampler::saveloops(std::ostream & os, const std::vector<NSPgeometry::XYZ> &crds,
		const std::array<double,IntrctBlck::ENESIZE> & energies) const {
	for(int i=0;i<IntrctBlck::ENESIZE;i++){
		os<<" "<<energies[i];
	}
	os<<std::endl;
	for(int i=0;i<loopatoms_.size();++i){
		os<<crds[loopatoms_[i]].toString()<<std::endl;
	}
}
std::vector<int> chooseloops(NSPdstl::RandomEngine<> &rng,int ntotal){
	if(ntotal<2) return std::vector<int>(1,0);
	double r=rng.randomreal();
	std::vector<int> res;
	if(r<0.5) {
		res.push_back(rng.intrng(0,ntotal-1)());
		return res;
	}
	int nsel=rng.intrng(2,ntotal)();
	while(res.size()<nsel){
		int ir=rng.intrng(0,ntotal-1)();
		bool chosen=false;
		for(auto i:res){
			if(ir==i) {
				chosen=true;
				break;
			}
		}
		if(!chosen) res.push_back(ir);
	}
	return res;
}
void LoopSampler::runMCCycles(){
	NSPdstl::RandomEngine<> rng;
	rng.init(mypara_.mcrandomseed);
	rng.setrealrng(0.0,1.0);
	std::ofstream outloops;
	outloops.open(mypara_.outloopfile);
	StopJudger mcjudger=mypara_.mkmcjudger();
	std::array<double,IntrctBlck::ENESIZE> energies;
	double eprev=100000000.0;
	std::vector<NSPgeometry::XYZ> currentcrds=imol_->recollectcrds_all();
	int nsaved=0;
	int ncycles=0;
	int nloops=looplengths_.size();
	while(true){
		this->rebuildloops(chooseloops(rng,nloops));
		this->optimizeloops();
		double enew=imol_->sum_energies(energies);
		bool accept=false;
		if(enew <= eprev) {
			accept=true;
		} else {
			double de=(enew-eprev)/mypara_.mckbt;
			if(de<4.0) if(rng.randomreal()<exp(-de)) accept= true;
		}
		if(accept){
			currentcrds=imol_->recollectcrds_all();
			eprev=enew;
			outloops<<nsaved++<<std::endl;
			saveloops(outloops,currentcrds,energies);
		} else {
			imol_->changecrds(currentcrds);
		}
		bool tostop=mcjudger(ncycles,enew);
		if(mypara_.verbose){
			std::cout <<"Minimum energy persisted for " <<mcjudger.eminpersist()<<" cycles"<< std::endl;
		}
		if(tostop) break;
		ncycles++;
	}
}
void NSPsampling::readpara_lpsampler(const std::string &parfile,
		LpSamplerPara &lpspara, NSPintrct::IntrctPara &ipara, SDRunPara &spara,
		std::string &jobname){
	NSPdataio::ControlFile cf;
	cf.readfile(parfile);
	lpspara=makelpsamplerparam(std::string(),cf);
	if(jobname=="auto"){
		jobname=lpspara.jobname;
	} else {
		lpspara.jobname=jobname;
	}
	ipara=makeintrctparam(jobname,cf);
	spara=makesdrunparam(jobname,cf);
	return;
}
LoopSampler NSPsampling::mkloopsampler(const LpSamplerPara &lpspara,
		const NSPintrct::IntrctPara &ipara,
		const SDRunPara &sdpara){
	NSPproteinrep::AAConformersInModel model,model1;
	typedef NSPintrct::MolModeler::SiteIdx SiteIdx;
	std::vector<std::pair<SiteIdx, SiteIdx>> findsflankingsites;
	std::vector<int> findlooplengths;
	std::pair<std::vector<std::pair<SiteIdx, SiteIdx>>, std::vector<std::pair<std::pair<int, int>, int>>> sitesandlength;
	model.readpdbfile(lpspara.pdbstart);
	model1.readpdbfile(lpspara.pdbstart);
	IntrctMol mol0(model);
	IntrctMol mol1(model1);
	LoopSampler res;
	res.mypara()=lpspara;
	res.mypara().sdprintsteps=sdpara.printsteps;
	for(auto & fp:lpspara.flankingpositions){
		int cid=mol0.mappdbkeyint()->chainNumber(fp.first.first);
//		std::cout <<"Chain " <<cid<<std::endl;
		int left=mol0.mappdbkeyint()->
				posiNumber(NSPproteinrep::PdbReader::reskey(fp.first.second),fp.first.first);
//		std::cout <<"left " <<left<<std::endl;
		int right=mol0.mappdbkeyint()->
				posiNumber(NSPproteinrep::PdbReader::reskey(fp.second.second),fp.first.first);
//		std::cout << "right " << right << std::endl;
		res.mypara().flankingsites.push_back({{cid,left},{cid,right}});
	}
	//for (auto& fp : lpspara.flankingpositions) {
	//	std::cout <<"left"<< fp.first.second <<"right"<< fp.second.second << std::endl;
	//}
	NSPdstl::RandomEngine<>::getinstance().init(res.mypara().mcrandomseed);
	if (lpspara.loopsearchmode == "short") {
		std::cout << "loopsearchmode =short " << std::endl;
		sitesandlength = res.findshortloops(mol0, res.mypara().flankingsites, res.mypara().looplengths, lpspara, ipara, sdpara);
	}
	else if (lpspara.loopsearchmode == "fixshortloop") {
		std::cout << "loopsearchmode =fixshortloop " << std::endl;
		sitesandlength = res.fixshortloops(mol0, res.mypara().flankingsites, res.mypara().looplengths, lpspara, ipara, sdpara);
	}
	else if (lpspara.loopsearchmode == "fixloop") {
		std::cout << "loopsearchmode =fixloop " << std::endl;
		sitesandlength = res.fixloops(mol0, res.mypara().flankingsites, res.mypara().looplengths, lpspara, ipara, sdpara);
	}
	else { std::cout << "loopsearchmode =auto " << std::endl;sitesandlength = res.findloops(mol0, res.mypara().flankingsites, res.mypara().looplengths, lpspara, ipara, sdpara); }
	
	std::vector<int> looplengths;
	for (auto& a : sitesandlength.second) {
		looplengths.push_back(a.second);
	}
	std::vector<std::pair<int, int>> lrshift;
	for (auto& a : sitesandlength.second) {
		lrshift.push_back(a.first);
	}
	res.setupmolfind(mol1, sitesandlength.first, looplengths);
	assert(sitesandlength.first.size() == res.mypara().flankingsites.size());
	std::cout << "findlength" << std::endl;
	for (auto a : looplengths) {
		std::cout<<a << std::endl;
	}
	std::cout << "res.mypara().looplengths" << std::endl;
	for (auto a : res.mypara().looplengths) {
		std::cout << a << std::endl;
	}
	std::vector<std::string> newloops;
	assert(lpspara.flankingpositions.size() == looplengths.size());
	for (int i = 0;i < looplengths.size();i++) {
		auto fp = lpspara.flankingpositions[i];
		int lshift = lrshift[i].first;
		int rshift = lrshift[i].second;
		//std::cout << "left" << fp.first.second << "right" << fp.second.second << std::endl;
		//std::cout << "chainid " << lpspara.flankingpositions[i].first.first << std::endl;
		std::string chainid;
		chainid.push_back(lpspara.flankingpositions[i].first.first);
		newloops.push_back(chainid +" ");//chainid
		newloops.push_back(std::to_string(std::stoi(fp.first.second) + lshift) + " ");
		newloops.push_back(std::to_string(std::stoi(fp.second.second) + rshift) + " ");
		newloops.push_back(std::to_string(looplengths[i]) + " ");
	}
	//for (int i = 0;i < looplengths.size();i++) {
	//	newloops.push_back("A ");
	//	newloops.push_back(std::to_string(sitesandlength.first[i].first.second + 1)+" " );
	//	newloops.push_back(std::to_string(sitesandlength.first[i].second.second + 1)+" ");
	//	newloops.push_back(std::to_string(sitesandlength.second[i]) + " ");
	//}
	std::ofstream ofs("newloop");
	for (auto s : newloops) ofs << s;
	ofs << std::endl;
	ofs.close();

	res.imol()->my_intrctpara=std::shared_ptr<IntrctPara>(new IntrctPara(mkrescaledparam(ipara)));
	res.setupsdrun(sdpara);
	return res;
}

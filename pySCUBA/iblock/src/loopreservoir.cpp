/*
 * loopreservior.cpp
 *
 *  Created on: 2019年12月21日
 *      Author: hyiu
 */
#include "iblock/loopreservoir.h"
#include "backbone/backbonebuilder.h"
using namespace NSPintrct;
using namespace NSPproteinrep;

std::shared_ptr<LoopReservoir::LoopConf> LoopReservoir::popconf(int length){
	if(nterminal_|| cterminal_){
		reservoir_[length]=std::vector<std::shared_ptr<LoopConf>>();
		posi_[length]=1;
	}
	if(nterminal_){
	reservoir_[length].push_back(std::shared_ptr<LoopConf>(new LoopConf(BackBoneBuilder::buildbackwardbackbone(length,anchor_right_,
				std::vector<std::pair<int,int>>(),
				std::vector<std::pair<int,int>>(),
				std::set<int>()))));
	} else if(cterminal_){
		reservoir_[length].push_back(std::shared_ptr<LoopConf>(new LoopConf(BackBoneBuilder::buildforwardbackbone(length,anchor_left_,
				std::vector<std::pair<int,int>>(),
				std::vector<std::pair<int,int>>(),
				std::set<int>()))));
	}
	if(reservoir_.find(length)!=reservoir_.end()){
		if(posi_.at(length)<reservoir_[length].size()) {
			++posi_.at(length);
			return reservoir_.at(length).at(posi_.at(length)-1);
		}else {
			reservoir_[length].clear();
		}
	}
	int ntry=0;
	while (reservoir_[length].empty() && ntry++<max_tryclosing_){
		reservoir_[length]=BackBoneBuilder::buildlinkers(length,anchor_left_,
				anchor_right_,
				std::vector<std::pair<int,int>>(),
				std::vector<std::pair<int,int>>(),
				std::set<int>());
	}
	if(reservoir_[length].empty()) return std::shared_ptr<LoopReservoir::LoopConf>
	                  (new LoopReservoir::LoopConf());
	posi_[length]=1;
	return reservoir_.at(length).at(posi_[length]-1);
}
std::shared_ptr<LoopReservoir::LoopConf> LoopReservoir::popconf(int length, bool helixinloop){
	if(nterminal_|| cterminal_){
		reservoir_[length]=std::vector<std::shared_ptr<LoopConf>>();
		posi_[length]=1;
	}
	if(nterminal_){
	reservoir_[length].push_back(std::shared_ptr<LoopConf>(new LoopConf(BackBoneBuilder::buildbackwardbackbone(length,anchor_right_,
				std::vector<std::pair<int,int>>(),
				std::vector<std::pair<int,int>>(),
				std::set<int>()))));
	} else if(cterminal_){
		reservoir_[length].push_back(std::shared_ptr<LoopConf>(new LoopConf(BackBoneBuilder::buildforwardbackbone(length,anchor_left_,
				std::vector<std::pair<int,int>>(),
				std::vector<std::pair<int,int>>(),
				std::set<int>()))));
	}
	if(reservoir_.find(length)!=reservoir_.end()){
		//test
		if (!(reservoir_[length].empty())) {
			if (posi_.at(length) < reservoir_[length].size()) {
				++posi_.at(length);
				return reservoir_.at(length).at(posi_.at(length) - 1);
			}
			else {
				reservoir_[length].clear();
			}
		}
		

	}
	else
	{
		reservoir_[length] = std::vector<std::shared_ptr<LoopConf>>();
	}
	int ntry=0;
	//test
	if (reservoir_.find(length) == reservoir_.end()) { std::cout << "length not in reservoir_ " << std::endl; }
	while (reservoir_[length].empty() && ntry++<max_tryclosing_){
		if (helixinloop)
		{
			std::vector<std::pair<int, int>> helixregion;
			int i = (rand()%(length));
			int j = (rand()%(length-i))+i;
			helixregion.push_back(std::pair<int, int>(i,j));
			reservoir_[length]=BackBoneBuilder::buildlinkers(length,anchor_left_,
					anchor_right_,
					helixregion,
					std::vector<std::pair<int,int>>(),
					std::set<int>());
		}
		else
			reservoir_[length]=BackBoneBuilder::buildlinkers(length,anchor_left_,
					anchor_right_,
					std::vector<std::pair<int,int>>(),
					std::vector<std::pair<int,int>>(),
					std::set<int>());
	}
	if(reservoir_[length].empty()) return std::shared_ptr<LoopReservoir::LoopConf>
	                  (new LoopReservoir::LoopConf());
	posi_[length]=1;
	return reservoir_.at(length).at(posi_[length]-1);
}

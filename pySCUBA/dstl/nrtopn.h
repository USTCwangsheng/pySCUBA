/*
 * nrtopn.h
 *
 *  Created on: 2020年10月4日
 *      Author: hyiu
 */

#ifndef DSTL_NRTOPN_H_
#define DSTL_NRTOPN_H_
#include "dstl/topn.h"
#include <vector>
#include <algorithm>
namespace NSPdstl{

/**
 * non-redundant topn storage
 */
template<typename Obj>
struct NRTopN{
	int N{0};
	double similaritycut{0};
	NRTopN(){;}
	NRTopN(int n, double scut):N(n),similaritycut(scut){;}
	std::vector<ScoredType<Obj>> stored;
	template<typename Similarity>
	int findindex (Obj & obj,double score,const Similarity &sim){
		int toreplace{-1};
		double similarity_min=1.e6;
		if(stored.size()<N) toreplace=stored.size();
		else if(score >stored.back().score)return -1;
		else toreplace= stored.size()-1;
		for(int i=stored.size()-1;i>=0;--i){
			double s=sim(obj,stored[i].object);
			if(s<similaritycut){
				if(score >stored[i].score){
					toreplace=-1;
					break;
				}
				if(s<similarity_min) {
						toreplace=i;
						similarity_min=s;
//						std::cout <<" To replace: " <<toreplace << " " <<s<<std::endl;
					}
				}
			}
		return toreplace;
	}
	template<typename Similarity>
    bool store(Obj &obj, double score, const Similarity &sim){
		 int idx=this->findindex(obj,score,sim);
//		 std::cout <<"NRTopN find idx: " <<idx<<std::endl;
		 if(idx>=0){
			 if(idx>=nstored()){
				 stored.push_back(ScoredType<Obj>(obj,score));
			 } else {
				 stored[idx]=ScoredType<Obj>(obj,score);
			 }
			 if(stored.size()>1) std::sort(stored.begin(),stored.end());
//			std::cout <<"config stored: " <<nstored()<<std::endl;
			 return true;
		 }
		 return false;
	}
	int nstored() const {return stored.size();}
	const Obj &getstoredobj(int i)  const {return stored[i].object;}
	ScoredType<Obj> &getstored(int i)  {return stored[i];}
	bool trystore(double score) const {
		if (stored.size()<N) return true;
		return (score <stored.back().score);
	}
};
}



#endif /* DSTL_NRTOPN_H_ */

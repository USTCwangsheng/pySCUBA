/*
 * vec2d.h
 *
 *  Created on: 2019年12月1日
 *      Author: hyiu
 */

#ifndef DSTL_VEC2D_H_
#define DSTL_VEC2D_H_
#include <vector>
#include <cassert>
namespace NSPdstl{
typedef std::pair<int,int> Idx2D; //2D index,
template <typename T>
struct Vec2D{
	struct Idx2DSuit{
		typedef std::vector<int> IdxOffSet;
		IdxOffSet offsets;
		int totalsize{0};
		Idx2DSuit(){;}
		void update(const std::vector<std::vector<T>> &vv){
	    	int i=0;
	    	for(auto &v:vv) {offsets.push_back(i); i=i+v.size();}
	    	totalsize=i;
	    }
	    int idx1d(const Idx2D &id2)const {return offsets.at(id2.first)+id2.second;}

	    Idx2D idx2d(int id1) const{
	    	assert(id1<totalsize);
	    	for(int idx=0;idx<offsets.size()-1;++idx){
	    		if(id1<offsets[idx+1]){
	    			return std::make_pair(idx,id1-offsets[idx]);
	    		}
	    	}
	    	return std::make_pair(-1,-1);
	    }
	    T & getelement(std::vector<std::vector<T>> &vv, int id1) const{
	    	Idx2D id2=idx2d(id1);
	    	return vv[id2.first][id2.second];
	    }
	    const T & getelement(const std::vector<std::vector<T>> &vv, int id1) const{
	     	Idx2D id2=idx2d(id1);
	     	return vv.at(id2.first).at(id2.second);
	     }
	};
	std::vector<std::vector<T>> D_,backupD_;
	Idx2DSuit idxsuit;
	void backupdata() { backupD_ = D_; }
	void dropbackdata() { D_ = backupD_; }
	void cleanbackdata() { backupD_.clear(); }
	void updateidxsuit(){idxsuit=Idx2DSuit(D_);}
	Vec2D(const std::vector<std::vector<T>> &data) {
		D_=data;
		updateidxsuit();
	}
	Vec2D(){;}
	void push_backvec(const std::vector<T> & v){ D_.push_back(v);}
	void push_backele(int idxv,const T &o){D_[idxv].push_back(o);}
	std::vector<T> & vec(int i){return D_[i];}
	const std::vector<T> vec(int i) const {return D_.at(i);}
	const std::vector<T> &backvec() const {return D_.back();}
	std::vector<T> &backvec(){return D_.back();}
	const T & operator[](int i) const {return idxsuit.getelement(D_,i);}
    T & operator[](int i) {return idxsuit.getelement(D_,i);}
	const T & operator[](Idx2D i) const {return D_.at(i.first).at(i.second);}
    T & operator[](Idx2D i) {return D_.at(i.first).at(i.second);}
	const T & operator ()(int i,int j) const {return D_.at(i).at(j);}
	T &operator()(int i,int j) {return D_.at(i).at(j);}
	int sized1() const{return D_.size();}
	int sized2(int i) const{assert(i<D_.size());return D_.at(i).size();}
	int sizetotal() const{ return idxsuit.totalsize;}
	Idx2D idx2d(int id1) const{return idxsuit.idx2d(id1);}
	int idx1d(Idx2D id2) const{return idxsuit.idx1d(id2);}
	int idx1d(int i,int j) const{return idxsuit.idx1d(std::make_pair(i,j));}
	int offset(int i1) const {return idxsuit.offsets.at(i1);}
};
}

#endif /* DSTL_VEC2D_H_ */

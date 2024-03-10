/*
 * internalcrds.cpp
 *
 *  Created on: 2020年7月21日
 *      Author: hyiu
 */
#include "iblock/internalcrds.h"
#include "geometry/calculators.h"

using namespace NSPintrct;
using namespace NSPgeometry;

InternalCrds::ConnectionTable InternalCrds::getconnectiontable(
		const std::vector<XYZ> &crds, double rcut) {
	ConnectionTable ct(crds.size(), std::vector<int>());
	for (int i = 0; i < crds.size() - 1; ++i) {
		for (int j = i + 1; j < crds.size(); ++j) {
			double r = distance(crds.at(i), crds.at(j));
			if (r <= rcut) {
				ct[i].push_back(j);
				ct[j].push_back(i);
			}
		}
	}
	return ct;
}
std::vector<std::pair<int, int>> InternalCrds::getbonds(
		const InternalCrds::ConnectionTable &ct) {
	std::vector<std::pair<int, int>> bonds;
	for (int i = 0; i < ct.size(); ++i) {
		for (int j : ct.at(i)) {
			if (j > i)
				bonds.push_back(std::make_pair(i, j));
		}
	}
	return bonds;
}
std::vector<std::vector<int>> InternalCrds::getangles(
		const InternalCrds::ConnectionTable &ct) {
	std::vector<std::vector<int>> angles;
	for (int j = 0; j < ct.size(); ++j) {
		std::vector<int> a(3, -1);
		a[1] = j;
		for (int m = 0; m < ct.at(j).size() - 1; ++m) {
			int i = ct.at(j).at(m);
			for (int n = m + 1; n < ct.at(j).size(); ++n) {
				int k = ct.at(j).at(n);
				if (i < k) {
					a[0] = i;
					a[2] = k;
				} else if (k < i) {
					a[0] = k;
					a[2] = i;
				}
				angles.push_back(a);
			}
		}
	}
	return angles;
}

std::vector<std::vector<int>> InternalCrds::getimpdihs(
		const InternalCrds::ConnectionTable &ct) {
	std::vector<std::vector<int>> ips;
	for (int i = 0; i < ct.size(); ++i) {
		if (ct.at(i).size() < 3)
			continue;
		std::vector<int> a(4, -1);
		a[0] = i;
		a[1] = ct.at(i).at(0);
		a[2] = ct.at(i).at(1);
		for (int m = 2; m < ct.at(i).size(); ++m) {
			a[3] = ct.at(i).at(m);
			ips.push_back(a);
		}
	}
	return ips;
}

std::vector<std::vector<int>> InternalCrds::gettorsions(
		const InternalCrds::ConnectionTable &ct,
		const std::vector<std::pair<int, int>> & bonds) {
	std::vector<std::vector<int>> torsions;
	for(auto & bd:bonds){
		int j=bd.first;
		int k=bd.second;
		if(ct.at(j).size()<2) continue;
		if(ct.at(k).size()<2) continue;
		std::vector<int> t(4,0);
		t[1]=j;
		t[2]=k;
		for(int m=0;m<ct.at(j).size();++m){
			int i=ct.at(j).at(m);
			if(i==k) continue;
			t[0]=i;
			for(int n=0;n<ct.at(k).size();++n){
				int l=ct.at(k).at(n);
				if(l==j) continue;
				if(l==i) continue;
				t[3]=l;
				torsions.push_back(t);
			}
		}
	}
	return torsions;
}

/*
 * internalcrds.h
 *
 *  Created on: 2020年7月21日
 *      Author: hyiu
 */

#ifndef IBLOCK_INTERNALCRDS_H_
#define IBLOCK_INTERNALCRDS_H_
#include "geometry/xyz.h"
namespace NSPintrct{
class InternalCrds {
public:
	typedef std::vector<std::vector<int>> ConnectionTable;
	InternalCrds(const std::vector<NSPgeometry::XYZ> &crds,
			double rcut = 1.95) {
		auto ct = getconnectiontable(crds, rcut);
		init(ct);
	}
	InternalCrds(const ConnectionTable &ct) {
		init(ct);
	}
	void init(const ConnectionTable &ct) {
		bonds_ = getbonds(ct);
		angles_ = getangles(ct);
		impdihs_ = getimpdihs(ct);
		torsions_ = gettorsions(ct, bonds_);
	}
	const std::vector<std::pair<int, int>> & bonds() const {return bonds_;}
	const std::vector<std::vector<int>> & angles() const {return angles_;}
	const std::vector<std::vector<int>> & impdihs() const {return impdihs_;}
	const std::vector<std::vector<int>> & torsions() const {return torsions_;}

	static ConnectionTable getconnectiontable(
			const std::vector<NSPgeometry::XYZ> &crds, double rcut = 1.95);
	static std::vector<std::pair<int, int>> getbonds(const ConnectionTable &ct);
	static std::vector<std::vector<int>> getangles(const ConnectionTable &ct);
	static std::vector<std::vector<int>> getimpdihs(const ConnectionTable &ct);
	static std::vector<std::vector<int>> gettorsions(const ConnectionTable &ct,
			const std::vector<std::pair<int, int>> &bonds);
	static bool connected(const ConnectionTable &ct,int i,int j){
		for(auto & n:ct.at(i)){
			if( n==j) return true;
		}
		return false;
	}
private:
	std::vector<std::pair<int, int>> bonds_;
	std::vector<std::vector<int>> angles_;
	std::vector<std::vector<int>> impdihs_;
	std::vector<std::vector<int>> torsions_;
};
}




#endif /* IBLOCK_INTERNALCRDS_H_ */

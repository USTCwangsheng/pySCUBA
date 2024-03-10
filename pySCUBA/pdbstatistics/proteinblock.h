/*
 * proteinblock.h
 *
 *  Created on: 2017年6月13日
 *      Author: hyliu
 */

#ifndef PDBSTATISTICS_PROTEINBLOCK_H_
#define PDBSTATISTICS_PROTEINBLOCK_H_

#include <backbone/backbonesite.h>
#include <map>

namespace NSPpdbstatistics {
class ProteinBlock{
public:
	typedef std::map<std::string,std::vector<double>> ExtPBs;
	static std::map<char,std::vector<double>> stdpbs;
	template<typename IT>
	static std::vector<double> getpbtorsions(IT &iter,int length=5){
		std::vector<double> torsions;
		int lenh=length/2;
		torsions.push_back((iter-lenh)->psi());
		for(int i=-lenh+1;i<lenh;++i){
			torsions.push_back((iter+i)->phi());
			torsions.push_back((iter+i)->psi());
		}
		torsions.push_back((iter+lenh)->phi());
		return torsions;
	}
	static char pbtype(const std::vector<double> &dihs, double *dev=nullptr);
	static std::string extpbtype(const ExtPBs & extpbs,const std::vector<double> &dihs, double *dev=nullptr);
	static char pbtype(const std::vector<NSPproteinrep::BackBoneSite> &sites, int posi, double *dev=nullptr);
	static char pbtype(const std::vector<NSPproteinrep::BackBoneSite> &sites, bool keepss,int posi, double *dev=nullptr);
	static std::string extpbtype(const ExtPBs & extpbs,const std::vector<NSPproteinrep::BackBoneSite> &sites, int posi, double *dev=nullptr);
	static std::vector<double> deviations(char pbtype,const std::vector<double> &dihs);
	static std::vector<double> deviations(char pbtype,const std::vector<NSPproteinrep::BackBoneSite> &sites, int posi);
	static std::vector<double> extractdihvec(const std::vector<NSPproteinrep::BackBoneSite> &sites,int posi);
	static void setpbtypes(std::vector<NSPproteinrep::BackBoneSite> &chain,
			int posistart=-1,int posiend=-1);
	static void setpbtypes(std::vector<NSPproteinrep::BackBoneSite> &chain, bool keepss,
				int posistart=-1,int posiend=-1);
private:
	static double dihdiff(double d1,double d2){
		double diff=d1-d2;
		while(diff > 180.0) diff -=360.0;
		while (diff <-180.0) diff+=360.0;
		return diff;
	}
};
}



#endif /* PDBSTATISTICS_PROTEINBLOCK_H_ */

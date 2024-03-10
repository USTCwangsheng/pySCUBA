/*
 * nn_scorepep.h
 *
 *  Created on: 2017年8月9日
 *      Author: hyliu
 */

#ifndef PDBSTATISTICS_NN_PEPSCORER_H_
#define PDBSTATISTICS_NN_PEPSCORER_H_
#include "dstl/nnregressionmodel.h"
#include "dataio/datapaths.h"
#include <map>
namespace NSPpdbstatistics {
std::vector<double> triangencoder(const std::vector<double> &torsions,
		std::vector<std::vector<double>> *dcdt=nullptr);
class NN_PepScorer {
public:
	enum Type {
		COIL, HELIX, STRAND
	};
	template<Type TYPE>
	static NN_PepScorer & scorer() {
		static NN_PepScorer scorer;
		static std::map<Type, std::string> filenames { { COIL,
				"coilscorer.dat" }, { HELIX, "helixscorer.dat" }, { STRAND,
				"strandscorer.dat" } };
		if (!scorer.initialized_)
			scorer.setup(filenames.at(TYPE));
		return scorer;
	}
	static double energy(char pbtype,const std::vector<double> &pbtorsions,
			std::vector<double> *dedt=nullptr){
		double ene{0.0};
		if(pbtype=='m' || pbtype=='H'){
//			ene=scorer<HELIX>().nnestimator_(encoder<HELIX>(pbtorsions));
			return 0.0;
		} else if(pbtype=='d'||pbtype=='E'){
//			ene=scorer<STRAND>().nnestimator_(encoder<STRAND>(pbtorsions));
			return 0.0;
		} else {
			NN_PepScorer & sce=NN_PepScorer::scorer<COIL>();
			double esum=0.0;
			if(dedt){
				dedt->assign(pbtorsions.size(),0.0);
			}
			for (int i=0;i< sce.nnestimators_.size();++i){
				if(!dedt) {
					esum +=(sce.nnestimators_[i])(encoder(pbtorsions,
							sce.encoder_intervals_[i]));
				} else {
					std::vector<std::vector<double>> dcdt;
					std::vector<double> dedc;
					esum +=(sce.nnestimators_[i])(encoder(pbtorsions,
							sce.encoder_intervals_[i],&dcdt),dedc);
					int nbit=(int)(360.0001/sce.encoder_intervals_[i]);
					for(int l=0; l<pbtorsions.size();++l){
						for(int m=0;m<nbit;++m){
							(*dedt)[l] += dedc[l*nbit+m]*dcdt[l*nbit+m][l];
						}
					}
				}
			}
			ene=esum/(double)(sce.nnestimators_.size());
			if(dedt){
				for(auto &d:*dedt) d/=-(double)(sce.nnestimators_.size());
			}
		}
		return -ene;
	}
	static std::vector<double> encoder(const std::vector<double> &pbtorsions, double
			intv,std::vector<std::vector<double>> *dcdt=nullptr){
//		static std::map<Type,double> intervals{{COIL,40.0},{HELIX,40.0},{STRAND,40.0}};
//		static double intv=intervals[TYPE];
		int nbit=(int)(360.00001/intv);
		int nbit_total=pbtorsions.size()*nbit;
		std::vector<double> code(nbit_total,0.0);
		int start=0;
		if(dcdt){
			dcdt->assign(nbit_total,std::vector<double>(pbtorsions.size(),0.0));
		}
		int  tidx=0;
		for (auto t:pbtorsions){
			if(t<0) t+=360.0;
			else if(t>360.0) t -=360.0;
			int idx=(int) (t/intv);
		    if(idx == nbit) idx=nbit-1;
		    double x=(t-intv*(double)idx)/intv;
		    code[start+idx]=1-x;
		    if(idx ==nbit-1) code[start]=x;
		    else code[start+idx+1]=x;
		    if(dcdt){
		    	(*dcdt)[start+idx][tidx]=-1.0/intv;
		    	if(idx==nbit-1)(*dcdt)[start][tidx]=1.0/intv;
		    	else (*dcdt)[start+idx+1][tidx]=1.0/intv;
		    }
		    start +=nbit;
		    ++tidx;
		}
		return code;
	}
private:
	std::vector<NSPdstl::L3Estimator> nnestimators_;
	std::vector<double> encoder_intervals_;
	bool initialized_ { false };
	void setup(const std::string &filename) {
//		std::string fn = NSPdataio::datapath() + filename;
		std::string fn=NSPdataio::datafilename(filename);
		std::ifstream ifs;
		ifs.open(fn.c_str());
		char buffer[120];
		ifs.getline(buffer,120);
		int nmodels=std::stoi(std::string(buffer));
		encoder_intervals_.clear();
		for (int i=0;i<nmodels;++i){
			ifs.getline(buffer,120);
			double intv=std::stod(std::string(buffer));
			encoder_intervals_.push_back(intv);
		}
		nnestimators_.clear();
		for(int i=0;i<nmodels;++i) {
			nnestimators_.push_back(NSPdstl::readl3estimator(ifs));
		}
		initialized_=true;
	}
};
class NN_SSClassifier {
public:
	static const std::map<char,int> SSIDX;
//	static const std::map<char,int> ssidx{{'H',0},{'E',1},{'C',2}};
	static NN_SSClassifier & getinstance() {
		static NN_SSClassifier classifier;
		static std::string filenames {"ssclassifier.dat" };
		if (!classifier.initialized_)
			classifier.setup(filenames);
		return classifier;
	}
	std::vector<double> probabilities(const std::vector<double> &pbtorsions,
			std::vector<std::vector<double>> *dp3dt);
	double probability(char sstype,const std::vector<double> &pbtorsions,
			std::vector<double> *dpdt){
		std::vector<std::vector<double>> dp3dt;
		std::vector<double> p3=probabilities(pbtorsions,&dp3dt);
		*dpdt=dp3dt[SSIDX.at(sstype)];
		return p3[SSIDX.at(sstype)];
	}
private:
	NSPdstl::L3SoftMaxClassifier classifier_;
	bool initialized_{false};
	double encoder_intv_;
	void setup(const std::string & filename);
};
}

#endif /* PDBSTATISTICS_NN_PEPSCORER_H_ */

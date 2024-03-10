/*
 * nnterms.h
 *
 *  Created on: 2019年12月4日
 *      Author: hyiu
 */

#ifndef IBLOCK_NNTERMS_H_
#define IBLOCK_NNTERMS_H_

#include "iblock/phipsicodes.h"
#include "dstl/nnregressionmodel.h"
#include "dataio/datapaths.h"
//#define USE_CB
/**
 * codes of this name space handles molecular structure-dependent interactions
 */
namespace NSPintrct{

/**
 * base class of neuro-network energy terms
 */
class NNTerm{
public:
	NNTerm(){;}
	virtual NSPdstl::CombinedEstimator &getestimator()=0;
	virtual void setup(const std::vector<PhiPsiCodes> & phipsicodes, int posi){
		abort(); //this function should not be called;
	}

	/**
	 * evaluate energy and derivatives
	 * get called when derived NNterm objects yield a single value
	 */
	double outvalue(std::vector<DvDxi> *dvdxi);

	virtual void setup(const std::vector<double> &crd,
			const std::vector<BSInChain> &bsinchain1,
			const std::vector<PhiPsiCodes> & phipsicodes1, int posi1,
			const std::vector<BSInChain> &bsinchain2,
						const std::vector<PhiPsiCodes> & phipsicodes2, int posi2){
		setup(phipsicodes1,posi1);
	}

	virtual std::vector<double> outvalues (std::vector<std::vector<DvDxi>> *dvdxi){
			std::vector<double> res;
			dvdxi->push_back(std::vector<DvDxi>());
			res.push_back(outvalue(&(dvdxi->back())));
			return res;
	}
	virtual ~NNTerm(){;}
protected:
	std::vector<double>  codevec_;
	std::vector<std::vector<DvDxi>> dcdx_;
};

/**
 * the phipsi nnterm
 */
template <int DISTRTYPE=-1>
class PhiPsiNNTerm :public NNTerm{
public:
	enum DistrType {COIL,MIXCOIL,GLY,PRO,PREPRO,HELIX,STRAND,TRANSPRO};
	/**
	 * read the network from a file on first call
	 *
	 * return the network in first and later calls
	 */
	static NSPdstl::CombinedEstimator & gettheestimator(){
		static NSPdstl::CombinedEstimator lsestimator;
		static bool estimatorread { false };
		if (!estimatorread) {
#ifdef _OPENMP
#pragma omp critical(phpsinnterm_global)
			{
				if(!estimatorread){
#endif
			std::string filename;
			if(DISTRTYPE==MIXCOIL){
//				filename = NSPdataio::datapath() + "mixcoilphipsi_nnmodel.dat";
				filename = NSPdataio::datafilename("mixcoilphipsi_nnmodel.dat");
			}
			std::ifstream ifs;
			ifs.open(filename.c_str());
			lsestimator.setup(ifs);
			estimatorread = true;
#ifdef _OPENMP
#pragma omp flush
				}
			}
#endif
		}
		return lsestimator;
	}

	/**
	 * wrapper for the virtual method call
	 */
	NSPdstl::CombinedEstimator &getestimator(){
		return gettheestimator();
	}

	/**
	 * estimate energy and forces for a position
	 */
	double phipsiene(const std::vector<PhiPsiCodes> & phipsicodes, int posi,
			std::vector<DvDxi> *dedx){
		codevec_.clear();
		dcdx_.clear();
		for (int i=0;i<6;++i){
					codevec_.push_back(phipsicodes[posi].phicodes[i]);
					dcdx_.push_back(phipsicodes[posi].dphicodesdx[i]);
			}
		for (int i=0;i<6;++i){
					codevec_.push_back(phipsicodes[posi].psicodes[i]);
					dcdx_.push_back(phipsicodes[posi].dpsicodesdx[i]);
		}
		return outvalue(dedx);
	}
private:
};

/**
 * local structure energy term
 */
class LSNNTerm :public NNTerm{
public:
	static NSPdstl::CombinedEstimator & gettheestimator();
	NSPdstl::CombinedEstimator &getestimator(){
		return gettheestimator();
	}
	void setup(const std::vector<PhiPsiCodes> & phipsicodes, int posi);
private:
};
/**
 * deprecated codes
 */
class CBData{
public:
	CBData(){;}
	CBData(std::vector<NSPgeometry::XYZ> &crdncac);
	const NSPgeometry::XYZ &cbcrd() const{ return cbcrd_;}
	const std::vector<std::vector<NSPgeometry::XYZ>> &drcb() const {return drcb_;}
	std::vector<NSPgeometry::XYZ> distributederiv(const NSPgeometry::XYZ derivrcb) const;
private:
	NSPgeometry::XYZ cbcrd_;
	std::vector<std::vector<NSPgeometry::XYZ>> drcb_;
};

/**
 * mainchain-mainchain pair energy calculator
 */
class SitePairNNTerm :public NNTerm{
public:
	static NSPdstl::CombinedEstimator & gettheestimator();
	NSPdstl::CombinedEstimator &getestimator(){
		return gettheestimator();
	}
	virtual void setup(const std::vector<double> &crd,
			const std::vector<BSInChain> &bsinchain1,
			const std::vector<PhiPsiCodes> & phipsicodes1, int posi1,
			const std::vector<BSInChain> &bsinchain2,
			const std::vector<PhiPsiCodes> & phipsicodes2, int posi2);
	/**
	 * deprecated codes
	 */
	void setupwithcbdata(const std::vector<double> &crd,
				const std::vector<BSInChain> &bsinchain1, const CBData & cbdata1,
				const std::vector<PhiPsiCodes> & phipsicodes1, int posi1,
				const std::vector<BSInChain> &bsinchain2, const CBData &cbdata2,
				const std::vector<PhiPsiCodes> & phipsicodes2, int posi2);
	double rca() const {return rca_;}
	const std::vector<DvDxi> & drcadx() const {return drcadx_;}
private:
	void add_dmcodevec(const std::vector<double> &crd,
			const BSInChain &bs1,const BSInChain &bs2);
	void add_dmcodevec(const std::vector<double> &crd,
			const BSInChain &bs1,const CBData & cbdata1,
			const BSInChain &bs2,const CBData &cbdata2);
	std::vector<double> encode_r(double r,std::vector<double> *dcdr);
	double rca_;
	std::vector<DvDxi> drcadx_;
};
/**
 * secondary structure terms, used for optional secondary structure type restraints
 */
class NN_SSTerm {
public:
	enum SSIDX{HELIX,STRAND,COIL,TERMINUS};
	static int sstype(std::vector<double> probabilities){
		if(probabilities.empty()) return TERMINUS;
		int sidx=0;
		for(int i=1;i<3;++i){
			if(probabilities[i]>probabilities[sidx])sidx=i;
		}
		return sidx;
	}
	static NSPdstl::L3SoftMaxClassifier & getclassifier();
	std::vector<double> probabilities(const std::vector<PhiPsiCodes> & phipsicodes,int posi,
			std::vector<std::vector<DvDxi>> *dp3dx);
	double probability(int ssidx,const std::vector<PhiPsiCodes> &phipsicodes,int posi,
			std::vector<DvDxi> *dpdx){
		std::vector<std::vector<DvDxi>> dp3dx;
		std::vector<double> p3=probabilities(phipsicodes,posi,&dp3dx);
		*dpdx=dp3dx[ssidx];
		return p3[ssidx];
	}
private:
	NSPdstl::L3SoftMaxClassifier classifier_;
	void setup(const std::string & filename);
};
/**
 * backbone-dependent rotamer conformational energy term
 */
class ConformerCode;
class NN_KaiTerm:public NNTerm{
public:
	static NSPdstl::CombinedEstimator &gettheestimator(const std::string &restype);
	NSPdstl::CombinedEstimator &getestimator(){
		return gettheestimator(resname_);
	}
	NN_KaiTerm(){;}
	NN_KaiTerm(const ConformerCode &conformercode);
private:
	std::string resname_;
};
/**
 * backbone-independent side chain conformer terms
 *
 * applied to terminus residues
 */
class NN_KaiTerm_T:public NNTerm{
public:
	static NSPdstl::CombinedEstimator &gettheestimator(const std::string &restype);
	NSPdstl::CombinedEstimator &getestimator(){
		return gettheestimator(resname_);
	}
	NN_KaiTerm_T(){;}
	NN_KaiTerm_T(const ConformerCode &conformercode);
private:std::string resname_;
};
/**
 * local mainchain-mainchain hydrogen bond term.
 * only apply between positions of short sequence separations
 */
class LocalBbHBNNTerm : public NNTerm
{
public:
    static NSPdstl::CombinedEstimator & gettheestimator();
    NSPdstl::CombinedEstimator &getestimator(){
        return gettheestimator();
    }
    virtual void setup(const std::vector<std::vector<double>> &dmatrix,
                       const std::vector<int> &atomids);
    virtual void setup(const std::vector<NSPgeometry::XYZ> &crds,
                        const std::vector<int> &atomids);
    virtual void setup(const std::vector<double> &crd,
            const std::vector<BSInChain> &bsinchain1, int posi1,
            const std::vector<BSInChain> &bsinchain2, int posi2);
    double rnomin() const {return rnomin_;}
    const std::vector<DvDxi> & drcdx() const {return drcdx_;}
private:
    void add_dmcodevec(const std::vector<std::vector<double>> &dmatrix,
                       const std::vector<int> &atomids);
    void add_dmcodevec(const std::vector<NSPgeometry::XYZ> &crd,
                       const std::vector<int> &atomids);
    void add_dmcodevec(const std::vector<double> &crd,
            const BSInChain &bs1,const BSInChain &bs2);
    std::vector<double> encode_r(double r,std::vector<double> *dcdr);
    double rnomin_;
    std::vector<DvDxi> drcdx_;
};

class LocalBbHBGeoNNTerm : public NNTerm
{
public:
    static NSPdstl::CombinedEstimator & gettheestimator();
    NSPdstl::CombinedEstimator &getestimator(){
        return gettheestimator();
    }
    double outvalue(std::vector<DvDxi> *dvdxi) {
        double result = NNTerm::outvalue(dvdxi);
        for (auto & d : *dvdxi) {
            d.first = atomids_[d.first];
            d.second = srno_ * d.second;
        }
        for (auto & d : dsrdx_) {
            dvdxi->push_back(std::make_pair(atomids_[d.first], result * d.second));
        }
        return result * srno_;
    }
    virtual void setup(const std::vector<NSPgeometry::XYZ> &crds,
                       const std::vector<int> &atomids) {
        std::vector<NSPgeometry::XYZ> internalcrds;
        for (int id : atomids) {
            internalcrds.push_back(crds[id]);
        }
        setup(internalcrds);
        atomids_ = atomids;
    }
    virtual void setup(const std::vector<double> &crd,
            const std::vector<BSInChain> &bsinchain1, int posi1,
            const std::vector<BSInChain> &bsinchain2, int posi2);
    double rnomin() const {return rnomin_;}
    const std::vector<DvDxi> & drcdx() const {return drcdx_;}
private:
    std::vector<double> encode_r(double r, std::vector<double> *dcdr);
    virtual void setup(const std::vector<NSPgeometry::XYZ> &crds);
    int calc_geometry(const std::vector<NSPgeometry::XYZ> &coords,
                      std::vector<double> &geometry,
                      std::vector<std::vector<DvDxi>> &dgdxi);
    void make_codevec(const std::vector<double> &geometry,
                      std::vector<std::vector<DvDxi>> &dgdxi);
    double radian2degree(const double rad);
    std::vector<int> atomids_;
    double srno_;
    std::vector<DvDxi> dsrdx_;
    double rnomin_;
    std::vector<DvDxi> drcdx_;
    inline double switch1d(double x0, double x1, double x, double *deriv) {
        *deriv = 0.0;
        if (x <= x0)
            return 1.0;
        else if (x >= x1)
            return 0.0;
        if (x <= x0 + 0.5 * (x1 - x0)) {
            double xs = (x - x0) / (x1 - x0);
            double fx = 1.0 - 2.0 * xs * xs;
            *deriv = -4.0 * xs/(x1-x0);
            return fx;
        } else {
            double xs = (x1 - x) / (x1 - x0);
            double fx = 2.0 * xs * xs;
            *deriv = -4.0 * xs/(x1-x0);
            return fx;
        }
    }
};

}





#endif /* IBLOCK_NNTERMS_H_ */

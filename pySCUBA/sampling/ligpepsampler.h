/*
 * ligpepsampler.h
 *
 *  Created on: 2020年10月23日
 *      Author: hyiu
 */

#ifndef LIGPEPSAMPLER_H_
#define LIGPEPSAMPLER_H_
#include "sampling/ligandpep.h"
#include "sampling/restraints.h"
#include "sampling/samplingutil.h"
#include "sampling/molsdrun.h"
namespace NSPsampling{
struct LigPSamplerPara{
		std::string jobname{""};
		std::string pdbstart{""};   //receptor PDB filename
	    std::string ifresidues;  //interface receptor residues
	    std::string saveconfigfile{""};
	    int ligpeplength; //length of designed ligand peptide
	    double ncontacts0;  //minimum number of ligand residues in contact with receptor interface residues
	    double kcontactres; //forcecontant for ligand-receptor contact restraints
		int randomseed{1357}; //<* random seed for randomengine of the MC cycles;
		int sdprintsteps{100}; //<*will be taken from the sdpara
		bool verbose{true};
		int maxsdsteps{0}; //<* maximum number of SD steps to run for quenching the loops in each cycle;
		double enedecay{0.99}; //<* per-step decay for estimating energy variations by time averaging;
		double sdenevarcut{10.0}; //<*Quenching SD stops when timer-averaged energy variations below this cutoff
		LigPSamplerPara(const std::vector<std::string> &controllines=std::vector<std::string>());
		StopJudger mksdjudger() const {
			return StopJudger(StopJudger::VAR,maxsdsteps,enedecay,sdenevarcut);
		}
};
/**
 * read LisPepSampler controlling paratemers from control file
 */
inline LigPSamplerPara makeligpsamplerparam(const std::string &jobname,const NSPdataio::ControlFile &cf,
		const std::string & controlname="LigPepSamplerPar"){
		std::vector<std::string> lines=cf.getcontrolines(controlname);
		if(!jobname.empty()) lines.push_back("JobName = "+jobname);
		return LigPSamplerPara(lines);
}
/**
 * read control file defining the sets of parameters needed by LigPepSampler
 */
void readpara_ligpepsampler(const std::string &parfile,
		LigPSamplerPara &lpspara, NSPintrct::IntrctPara &ipara, SDRunPara &sdpara,
		std::string &jobname);
class LigPepSampler{
public:
	LigPepSampler(){;}
	void setupmol(const NSPintrct::IntrctMol & receptor);
	/**
	 * Access the molecule hosting the loops
	 */
	std::shared_ptr<NSPintrct::IntrctMol> imol() const{return imol_;}
	void setupsdrun(const SDRunPara &sdpara);
	void rebuildligand(){
		this->ligandpep_->genligandpepcrd(imol_.get(),imol_->nchains()-1);
	}
	void buildandoptimize();
	void saveconfig(std::ofstream & ofs) const;
	LigPSamplerPara & mypara() {return mypara_;}
	const LigPSamplerPara & mypara() const {return mypara_;}
private:
	LigPSamplerPara mypara_;
	std::shared_ptr<NSPintrct::IntrctMol> imol_{nullptr}; //*< The IntrctMol comprised of receptor and loop
	std::shared_ptr<LigandPep> ligandpep_{nullptr};
	MolSDRun optimizer_; //*< for conformation optimization
};
/**
 * construct loop sampler based on the controlling parameters
 */
LigPepSampler mkligpepsampler(const LigPSamplerPara &lpspara,
		const NSPintrct::IntrctPara &ipara,
		const SDRunPara &sdpara);
}


#endif /* LIGPEPSAMPLER_H_ */

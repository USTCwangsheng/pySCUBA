/*
 * loopsampling.h
 *
 *  Created on: 2019年12月23日
 *      Author: hyiu
 */

#ifndef LOOPSAMPLER_H_
#define LOOPSAMPLER_H_
#include "iblock/molmodeler.h"
#include "sampling/molsdrun.h"
#include "sampling/samplingutil.h"
#include "iblock/blkselector.h"
#include "iblock/intrctblck.h"
#include <stdlib.h>
#include <vector>
#include <random>
#include <algorithm>
namespace NSPsampling{
/**
 * Input parameters for setup and run the Loop Sampler
 *
 * The program testloopsampler illustrate how to use this class
 */
struct LpSamplerPara{
	typedef NSPintrct::MolModeler::SiteIdx SiteIdx;
	typedef std::pair<char,std::string> PdbPosiIdx;
	std::vector<std::pair<PdbPosiIdx,PdbPosiIdx>> flankingpositions;
	int nloops{0}; //*< number of loops;
	std::vector<std::pair<SiteIdx,SiteIdx>> flankingsites; //*< positions of loops as indexed in the starting protein
	std::vector<int> looplengths;
	std::string jobname{""};
	std::string pdbstart{""};
	std::string outloopfile{""};
	int mcrandomseed{1357}; //<* random seed for randomengine of the MC cycles;
	int sdprintsteps{100}; //<*will be taken from the sdpara
	bool verbose{true};
	double mckbt{10.0}; //<* KBT for Metropolis MC cycles
	int maxcycles{0}; //<* maximum number of loop regenerating cycles;
	int maxsdsteps{0}; //<* maximum number of SD steps to run for quenching the loops in each cycle;
	double enedecay{0.99}; //<* per-step decay for estimating energy variations by time averaging;
	double sdenevarcut{10.0}; //<*Quenching SD stops when timer-averaged energy variations below this cutoff
	int MCcyclescut{20}; //*<the outer MC cycle stops if the lowest energy does not change in a given number of cycles.
	std::string pmax{"auto"}; // Max possibility of choosing init loops from LoopPool. auto = 0.95
	std::string unchangedcut{"auto"}; // If continuously no new loops, then P = PMax. auto = max(1000, 5*ceil(pow(1.7, LoopLength_Average)).
	std::string poolsizecut{"auto"}; // P(choosing old loop) = PoolSize/PSC*PMax. auto = 5*ceil(pow(1.7, ll_ave)
	int reconstructnum{0}; // Num of runs you'd like to sample.
	bool helixinloop{false}; // random region in the sampled loop using Helix_phipsi as init conf.
	int readloopnum{0}; // Num of non-redundant loops you'd like to fetch.
	std::string readrmsd{"auto"}; // could be auto, min, or (double)RMSD.
	std::string phipsirmsd{"auto"}; // could be auto or (double)RMSD.
	std::string loopsearchmode{ "auto" };// could be auto or short.
	int findlooptimes{50};//The number of times used to calculate the average energy of the loop
	LpSamplerPara(const std::vector<std::string> &controllines=std::vector<std::string>());
	StopJudger mkmcjudger() const {
		return StopJudger(StopJudger::MINPERSIST,maxcycles,MCcyclescut,0.0);
	}
	StopJudger mksdjudger() const {
		return StopJudger(StopJudger::VAR,maxsdsteps,enedecay,sdenevarcut);
	}

};
/**
 * each LoopPool contains many loopconfs with their loop_Ene.
 */
struct LoopPool
{
	// for which loop
	std::map<double, std::vector<NSPintrct::IntrctBlck>> loopconfs; // loop_Ene, loop_conf.
	double p_frompool = 0.0; // posibility to choose
	int unchanged = 0; // if unchanged > UNC_th (e.g. 200), then p_frompool = 1.0;
	double P_max = 1.0;
	int Unchanged_th = 200;
	int Size_th = 20;
};
/**
 * Wraps MolModler and MolSDRun to construct and optimize loops in peptide chains
 *
 * it handles multiple loops simultaneously. In each cycle, one or more
 * loops will be constructed using random torsion angles drawn from coil phi-psi distributions,
 * analytically closed, and then quenched to local minima of SCUBA energy function using SD.
 * The resulting loop conformations will be written to a loop trajectory file, together
 * with the SCUBA energy components. If in Monte Carlo mode, the reconstructed loops
 * will be accepted or rejected using a Metropolis criteria.
 */
class LoopSampler{
public:
	typedef NSPintrct::MolModeler::LoopConfMaker ConfMaker;
	typedef NSPintrct::MolModeler::SiteIdx SiteIdx;
	typedef std::pair<SiteIdx,SiteIdx>FlankingSite;
	/**
	 * Set up  work w.r.t the new molecule
	 *
	 * The host IntrctMol object is constructed based on the starting mol,
	 * loop sites, and new loop length. ConfMaker for each loop is constructed
	 *
	 * Returns false if failed to construct initial closed loops with the given
	 * loop flanking sites and loop lengths.
	 */
	bool setupmol(NSPintrct::IntrctMol &startmol,
			const std::vector<FlankingSite> & flankingsites,
		    const std::vector<int> & looplengths);
	bool setupmolfind(NSPintrct::IntrctMol& startmol,
		const std::vector<FlankingSite>& flankingsites,
		std::vector<int>& looplengths);
	/**
	 * set up work w.r.t. the quenching SD simulations to optimize the loops
	 */
	void setupsdrun(const SDRunPara & para){
		imol_->specifyactiveblks(inloopblks(reindexedflankingsites_),
				typename NSPintrct::BlkSelector::SelectedBlks());
		optimizer_.init(para,imol_);
	}

	NSPintrct::BlkSelector::SelectedBlks loopblks(const std::vector<std::pair<int, int>>& lranges);

	/**
	 * find best loop sites&lengths before setupmol
	 */
	std::pair<std::vector<std::pair<SiteIdx, SiteIdx>>, std::vector<std::pair<std::pair<int, int>, int>>> findloops(NSPintrct::IntrctMol& startmol,
		const std::vector<FlankingSite>& flankingsites,
		std::vector<int>& looplengths, const LpSamplerPara& lpspara,
		const NSPintrct::IntrctPara& ipara,
		const SDRunPara& sdpara);

	/**
	 * find short loop sites&lengths before setupmol
	 */
	std::pair<std::vector<std::pair<SiteIdx, SiteIdx>>, std::vector<std::pair<std::pair<int, int>, int>>> findshortloops(NSPintrct::IntrctMol& startmol,
		const std::vector<FlankingSite>& flankingsites,
		std::vector<int>& looplengths, const LpSamplerPara& lpspara,
		const NSPintrct::IntrctPara& ipara,
		const SDRunPara& sdpara);
	/**
	 * find fixed short loop lengths
	 */
	std::pair<std::vector<std::pair<SiteIdx, SiteIdx>>, std::vector<std::pair<std::pair<int, int>, int>>> fixshortloops(NSPintrct::IntrctMol& startmol,
		const std::vector<FlankingSite>& flankingsites,
		std::vector<int>& looplengths, const LpSamplerPara& lpspara,
		const NSPintrct::IntrctPara& ipara,
		const SDRunPara& sdpara);

	/**
	* find fixed  loop lengths
	*/
	std::pair<std::vector<std::pair<SiteIdx, SiteIdx>>, std::vector<std::pair<std::pair<int, int>, int>>> fixloops(NSPintrct::IntrctMol& startmol,
		const std::vector<FlankingSite>& flankingsites,
		std::vector<int>& looplengths, const LpSamplerPara& lpspara,
		const NSPintrct::IntrctPara& ipara,
		const SDRunPara& sdpara);

	/**
	 * save the quenched loop configurations with energy components
	 */
	void saveloops(std::ostream &os, const std::vector<NSPgeometry::XYZ> &crds,
			const std::array<double, NSPintrct::IntrctBlck::ENESIZE> &energies) const;

	/**
	 * rebuild loops indexed by lidx in imol_
	 *
	 */
	bool rebuildloops(const std::vector<int> &lidx);

	/**
	 * reuse loops from looppool, indexed by lidy in imol_
	 *
	 */
	bool reuseloops(const std::vector<int> &lidy);

	/**
	 * using SD simulation to optimize the loop
	 */
	void optimizeloops();

	std::shared_ptr<std::vector<NSPgeometry::XYZ>>
		           optimizeloops(double &ene,std::array<double,NSPintrct::IntrctBlck::ENESIZE> &energies){
		optimizeloops();
		ene=potenergy(energies);
		return loopatomcrds();
	}
	double potenergy(std::array <double,NSPintrct::IntrctBlck::ENESIZE> &energies) const {
		 return imol_->sum_energies(energies);
	}
	std::vector<double> loopenes(
			std::vector<std::array<double,NSPintrct::IntrctBlck::ENESIZE>> & energies) const;

	std::shared_ptr<std::vector<NSPgeometry::XYZ>>  loopatomcrds() const{
		 auto & currentcrd =imol_->recollectcrds_all();
		auto res=std::shared_ptr<std::vector<NSPgeometry::XYZ>>(new std::vector<NSPgeometry::XYZ>);
			for(int i=0;i<loopatoms_.size();++i){
				res->push_back(currentcrd.at(loopatoms_[i]));
			}
			return res;
	}
	std::shared_ptr<std::vector<NSPgeometry::XYZ>>
	           regenloops(double &ene,std::array<double,NSPintrct::IntrctBlck::ENESIZE> &energies){
		std::vector<int> lidx;
		for(int i=0;i<looplengths_.size();++i) lidx.push_back(i);
		rebuildloops(lidx);
		return optimizeloops(ene, energies);
	}
	/**
	 * random_gen OR use LoopPool
	 */
	std::shared_ptr<std::vector<NSPgeometry::XYZ>>
	           regenuseloops(double &ene,std::array<double,NSPintrct::IntrctBlck::ENESIZE> &energies){
		std::vector<int> lidx, lidy,box;
		std::random_device rd;
		std::mt19937 rng(rd());
		
		if (looppools_.size() == 0)
			looppools_.resize(looplengths_.size());

		int repalcenum = (rand()% 3) + 1;
		std::cout << "repalcenum: " << repalcenum<< std::endl;
		for (int i = 0;i < looplengths_.size();++i) {
			box.push_back(i);
		}
		std::shuffle(box.begin(), box.end(), rng);
		box.resize(repalcenum);
		for(auto a:box){ std::cout << a << ' '; }
		std::cout << std::endl;

		//for(int i=0;i<looplengths_.size();++i)
		//{
		//	if (rand()/(double)(RAND_MAX) >= looppools_[i].p_frompool)
		//		lidx.push_back(i);
		//	else
		//		lidy.push_back(i);
		//}
		for (auto i:box)
		{
			if (rand() / (double)(RAND_MAX) >= looppools_[i].p_frompool) {
				lidx.push_back(i);//std::cout << "rebuildloops: " << i << std::endl;
			}
			else {
				lidy.push_back(i);//std::cout << "looppools: " << i << std::endl;
			}
		}
		
		rebuildloops(lidx);
		reuseloops(lidy);
		return optimizeloops(ene, energies);
	}

	/**
	 * update looppools_ by modeler_ and reindexedflankingsites_
	 */
	void updatelooppools(){
		for (int i = 0; i < reindexedflankingsites_.size(); i++)
		{
			// confs.
			std::vector<NSPintrct::IntrctBlck> newloop =
					modeler_.outputloop(reindexedflankingsites_[i].first,reindexedflankingsites_[i].second);
			// ene.
			double etot_loop = 0.0;
			int cid=reindexedflankingsites_[i].first.first;
			int rstart=reindexedflankingsites_[i].first.second-1;
			if(rstart<0) rstart=0;
			int rend=reindexedflankingsites_[i].second.second+1;
			if(rend>=imol_->nresidues(cid)) rend=imol_->nresidues(cid)-1;
			for(int r=rstart; r<=rend; ++r){
				for(int t=0;t<NSPintrct::IntrctBlck::ENESIZE;++t){
					etot_loop +=imol_->getblck(NSPdstl::Idx2D(cid,r)).energies[t];
				}
			}
			// RMSD_th
			double RMSD_th = 0.05; // 0.03 + 0.02 * (newloop.size() - 4); // same with SavedLoops.rmsdcut()
			bool accept = false;
			if (etot_loop < 0.0)
				accept = true; // simple threshold for acception.
			if (accept)
			{
				bool direct_add = true;
				for (auto it = looppools_[i].loopconfs.begin(); it != looppools_[i].loopconfs.end();)
				{
					// see if redundant
					std::vector<NSPgeometry::XYZ> newcrds, oldcrds;
					for (int r = 0; r < newloop.size(); r++)
					{
						newcrds.insert(newcrds.end(), newloop[r].crds.begin(), newloop[r].crds.end());
						oldcrds.insert(oldcrds.end(), it->second[r].crds.begin(), it->second[r].crds.end());
					}
					assert(newcrds.size() == oldcrds.size());
					if (rmsd(newcrds, oldcrds) < RMSD_th && rmsd(newcrds, oldcrds) > 0)
					{
						if (it->first > etot_loop) // replace old_loop
						{
							looppools_[i].loopconfs.erase(it++);
							direct_add = false;
							looppools_[i].loopconfs.insert(std::make_pair(etot_loop, newloop));
							double p_new = looppools_[i].P_max *
									(double)looppools_[i].loopconfs.size() / looppools_[i].Size_th;
							if (p_new < looppools_[i].p_frompool
									&& looppools_[i].unchanged < looppools_[i].Unchanged_th)
								looppools_[i].p_frompool = p_new; // one new conf. replace many old confs.
							looppools_[i].unchanged = 0;
						}
						else  // discard new_loop
						{
							looppools_[i].unchanged++;
							direct_add = false;
							++it;
						}
					}
					else
						++it;
				} // check unredundant for each existing conf.
				if (direct_add)
				{
					looppools_[i].loopconfs.insert(std::make_pair(etot_loop, newloop));
					if (looppools_[i].unchanged >= looppools_[i].Unchanged_th &&
							looppools_[i].p_frompool < looppools_[i].P_max)
						looppools_[i].p_frompool = looppools_[i].P_max;
					if (looppools_[i].p_frompool < looppools_[i].P_max)
						looppools_[i].p_frompool = looppools_[i].P_max *
								(double)looppools_[i].loopconfs.size() / looppools_[i].Size_th;
					looppools_[i].unchanged = 0;
				}
			}
			// test_cyx
			std::cout << "LoopPool " << i << " with p_frompool " << looppools_[i].p_frompool << " & size "
					<< looppools_[i].loopconfs.size() << " & unchanged " << looppools_[i].unchanged << std::endl;
		} // each loop
	}

	/*
	 * SCUBALoopSampler: sampler read parloop from lpspara.
	 */
	void readparams(LpSamplerPara lps)
	{
		looppools_.resize(looplengths_.size());
		int ll_ave = 0;
		for (auto ll : looplengths_)
			ll_ave += ll;
		ll_ave = ceil(ll_ave/looplengths_.size()); // average loop lengths
		for (auto &lp : looppools_)
		{
			if (lps.pmax == "auto")
				lp.P_max = 0.95;
			else
				lp.P_max = std::stod(lps.pmax);
			if (lps.unchangedcut == "auto")
				lp.Unchanged_th = 1000 > 5*ceil(pow(1.7, ll_ave)) ? 1000 : 5*ceil(pow(1.7, ll_ave));
			else
				lp.Unchanged_th = std::stoi(lps.unchangedcut);
			if (lps.poolsizecut == "auto")
				lp.Size_th = 5*ceil(pow(1.7, ll_ave));
			else
				lp.Size_th = std::stoi(lps.poolsizecut);
		}
		helixinloop_ = lps.helixinloop;
	};

	/**
	 * use params for looppools_
	 */
	void useparams(std::string parloop)
	{
		looppools_.resize(looplengths_.size());
		int ll_ave = 0;
		for (auto ll : looplengths_)
			ll_ave += ll;
		ll_ave = ceil(ll_ave/looplengths_.size()); // average loop lengths
		std::ifstream ifpar(parloop.c_str()); // gnp.gnp
		if(!ifpar.good()) {
			std::cout << "parloop.params file failure" << std::endl;
			exit(1);
		}
		while(true) {
			std::string line;
			line.clear();
			getline(ifpar, line);
			if(!ifpar.good()) break;
			if (line.size()==0) continue;
			std::vector<std::string> words;
			std::stringstream input(line);
			std::string word;
			while(input>>word) words.push_back(word);
			if (words[0] == "P_max" && words.size() == 3)
			{
				double p_max;
				if (words[2] == "auto")
					p_max = 0.95;
				else
					p_max = std::stod(words[2]);
				for (auto &lp : looppools_)
					lp.P_max = p_max;
				std::cout << "P_max = " << p_max << std::endl;
			}
			if (words[0] == "Unchanged_th" && words.size() == 3)
			{
				int unchanged_th;
				if (words[2] == "auto")
					unchanged_th = 1000 > 5*ceil(pow(1.7, ll_ave)) ? 1000 : 5*ceil(pow(1.7, ll_ave)); //10 * ll_ave + 50;
				else
					unchanged_th = std::stoi(words[2]);
				for (auto &lp : looppools_)
					lp.Unchanged_th = unchanged_th;
				std::cout << "Unchanged_th = " << unchanged_th << std::endl;
			}
			if (words[0] == "Size_th" && words.size() == 3)
			{
				int size_th;
				if (words[2] == "auto")
					size_th = 5*ceil(pow(1.7, ll_ave)); //ceil(pow(1.5, ll_ave));
				else
					size_th = std::stoi(words[2]);
				for (auto &lp : looppools_)
					lp.Size_th = size_th;
				std::cout << "Size_th = " << size_th << std::endl;
			}
		}
		ifpar.close();
	}

	/**
	 * major MC cycle includes generating new loop configurations, quenching by SD and accept/reject
	 * the newly generated loop configuration
	 */
	void runMCCycles();
	/**
	 * determine IntrctBlcks forming the loops
	 */
	static NSPintrct::BlkSelector::SelectedBlks inloopblks(
			const std::vector<FlankingSite> &flankingsites);
	LpSamplerPara & mypara(){return mypara_;}
	const LpSamplerPara &mypara() const{return mypara_;}

	/**
	 * Access the molecule hosting the loops
	 */
	std::shared_ptr<NSPintrct::IntrctMol> imol() const{return imol_;}
	std::vector<FlankingSite> & reindexedflankingsites() {return reindexedflankingsites_;}
	const std::vector<FlankingSite> & reindexedflankingsites() const {return reindexedflankingsites_;}
private:
	std::shared_ptr<NSPintrct::IntrctMol> imol_{nullptr}; //*< The new IntrctMol that hosts the loops
	NSPintrct::MolModeler modeler_; //*<modeler targeting imol
	std::vector<ConfMaker> confmakers_; //*< initial loop conformation producer
	std::vector<FlankingSite> reindexedflankingsites_;
	std::vector<int> loopatoms_;
	std::vector<int> looplengths_;
	MolSDRun optimizer_; //*< for conformation optimization
	LpSamplerPara mypara_;
	std::vector<LoopPool> looppools_; // *< contain low_Ene loop confs. for each loop.
	bool helixinloop_{false}; // random region in the sampled loop using Helix_phipsi as init conf.
};
/**
 * read LoopSampler controlling paratemers from control file
 */
inline LpSamplerPara makelpsamplerparam(const std::string &jobname,const NSPdataio::ControlFile &cf,
		const std::string & controlname="LoopSamplerPar"){
		std::vector<std::string> lines=cf.getcontrolines(controlname);
		if(!jobname.empty()) lines.push_back("JobName = "+jobname);
		return LpSamplerPara(lines);
}

/**
 * read control file defining the sets of parameters needed by LoopSampler
 */
void readpara_lpsampler(const std::string &parfile,
		LpSamplerPara &lpspara, NSPintrct::IntrctPara &ipara, SDRunPara &spara,
		std::string &jobname);

/**
 * construct loop sampler based on the controlling parameters
 */
LoopSampler mkloopsampler(const LpSamplerPara &lpspara,
		const NSPintrct::IntrctPara &ipara,
		const SDRunPara &spara);
/**
 * determine whether chain of s1 is after chain of s2;
 * or if they are in the same chain, whether position of s1 if after position s2
 * To be used in determining the order of inserting loops.
 * If s1<s2, then inserting loop at s1 of any lengths will not affect the index of s2

inline bool operator<(const LoopSampler::FlankingSite &s1,
		const LoopSampler::FlankingSite &s2){
	if(s2.first.first<s1.first.first) return true;
	if(s2.first.second<s1.first.second) return true;
	return false;
}*/
/**
 * get new indices of the loop flankingsites after looplength change
 */
std::vector<LoopSampler::FlankingSite> reindexsites(const std::vector<LoopSampler::FlankingSite> & oldsites,
		const std::vector<int> & looplengths);
}
#endif /* LOOPSAMPLER_H_ */

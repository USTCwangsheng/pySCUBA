/*
 * pocketdocker.h
 *
 *  Created on: 2020年7月22日
 *      Author: hyiu
 */

#ifndef POCKETDOCKER_H_
#define POCKETDOCKER_H_
#include "proteinrep/pocketmodeler.h"
#include "iblock/intrctmol.h"
#include "geometry/rectgrid.h"
#include "geometry/sphere.h"
#define MACHTING_DISTANCE2MAX 4.0
namespace NSPsampling{
struct NbrFndrGrid{
	NSPgeometry::RectGrid grid;//<*grid to accelerate clash detection and so on
	typedef NSPgeometry::RectGrid::PntIdx PntIdx;
	std::map<PntIdx, std::vector<int>> gridneighbors;
	double cutoff2;
	void setgrid(const NSPgeometry::XYZ &center, const std::array<int,3> &sizes, double step);
	void setgridneighbors(const NSPintrct::IntrctMol & imol,double rcut2);
	std::vector<int> candidateneighbors(const NSPgeometry::XYZ &x) const {
		PntIdx xp;
		 for(int i=0;i<3;++i)
			 xp[i]=(x[i]-grid.p000[i])/grid.step;
		 auto entry=gridneighbors.find(xp);
		 if(entry!=gridneighbors.end()) return entry->second;
		 return std::vector<int>();
	}
};
struct BckBnStMatcher{
	static double sitedistance2(const NSPintrct::IntrctBlck &iblck,
				const NSPproteinrep::AAConformer &aaconf);
	NSPdstl::Idx2D matchsite(const NSPintrct::IntrctMol & imol,
			const NSPproteinrep::AAConformer &aaconf,double &score) const;
	std::set<NSPdstl::Idx2D> candidatesites; //<* candidate site for installing pocket residues
	double maxdist2{MACHTING_DISTANCE2MAX};
};
struct PocketMover{
		NSPgeometry::Sphere ligandcentersphere;
		double maxstep_trans{1.0};
		double maxstep_rotate{10}; //*< in degree
		double prigidmv{0.5};
		std::map<NSPproteinrep::PocketModeler::AtomID,NSPgeometry::XYZ> refcrd;
		std::vector<NSPgeometry::XYZ> crd_before;
		std::set<NSPproteinrep::PocketModeler::ComponentID> movencomponents;
		int pocketclashes;//<* current number of internal clashes within the isolated pocket
		void stepback(NSPproteinrep::PocketModeler &pm);
		bool randommv(NSPproteinrep::PocketModeler &pm);
		bool mvtosphere(NSPproteinrep::PocketModeler &pm);
		bool mvtofit(NSPproteinrep::PocketModeler &pm);
};
struct PckDckrParam{
	  enum POCKETLOCATIONMODE{CENTERXYZ, CENTERATOMS,REFATOMS};
	  std::string proteinpdb;
	  std::string pocketfile;
	  std::string ligandlocationmode;
	  NSPgeometry::XYZ ligandsphere_center;
	  std::vector<std::string> ligandcenteratoms;
	  double ligandsphereradius;
	  NSPgeometry::XYZ  gridcenter;
	  std::array<double,3> gridlengths;
	  double gridstep{3.0};
	  double maxtranslation{0.5};
	  double maxrotation{5.0};
};
PckDckrParam mk_pckdckrparam(const NSPdataio::ControlFile &cf);
class PocketDocker{
public:
	void setup(const PckDckrParam &param);
	bool micromove();
private:
	std::shared_ptr<NSPproteinrep::PocketModeler> pocketmodeler_; //<* isolated pocket and a sampler;
	std::shared_ptr<NSPintrct::IntrctMol> imol_; //<* receptor protein
	PocketMover pocketmover_;
	NbrFndrGrid nbrfndrgrid_;
	BckBnStMatcher bckbnstmatcher_;
	std::map<NSPproteinrep::IsolatedPocket::ComponentID,double> clash_scores_;
	             //<* current clash scores between pocket components and receptor protein atoms
	std::map<NSPproteinrep::IsolatedPocket::ComponentID,
	                std::pair<NSPdstl::Idx2D,double>> sitematchscores_;  //<* current matched pocket site and receptor site;
	double totalclashscore_;
	double totalmatchscore_;
	double score_clashes();
	double score_clashes(const NSPproteinrep::IsolatedPocket::ComponentID &comp);
	double score_sitematch();
	double score_sitematch(const NSPproteinrep::IsolatedPocket::ComponentID &comp);
};
}



#endif /* POCKETDOCKER_H_ */

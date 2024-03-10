/*
 * restraints.h
 *
 *  Created on: 2020年8月13日
 *      Author: hyiu
 */

#ifndef RESTRAINTS_H_
#define RESTRAINTS_H_
#include "geometry/calculators.h"
#include "iblock/ssrestraint.h"
#include "iblock/intrctmol.h"
#include "sampling/groupdistance.h"
namespace NSPsampling{
class ResTraint{
public:
	virtual  double energy(const std::vector<NSPgeometry::XYZ>& crd, std::vector<NSPgeometry::XYZ>* forces) const = 0;
	virtual ~ResTraint() { ; }
	double& ene() const { return ene_; }
//	int& switchable() const { return switchable_; }
	int switchable_ = 0;
	std::string restriantname_{"restriants"};
	std::string& restriantname() {return restriantname_;}
	int& switchable() { return switchable_; }
//	std::vector<int> switchkeys;

	double energy(const std::vector<NSPgeometry::XYZ>& crd, std::vector<NSPgeometry::XYZ>* forces, int resswitch) const {
		return ene_ = 0;
		}



	virtual void printinfo(std::ostream& os) const = 0;


protected:
	mutable double ene_{ 0.0 };
private:
	
};


class StructRestraint:public ResTraint{
public:
	enum {POSIMODE,TOTALMODE};
	StructRestraint(){;}
//	StructRestraint(const std::vector<NSPgeometry::XYZ> &crd, const std::vector<double> &weights,double kres,
//			int mode=POSIMODE,double rmsdref=0.0);
	StructRestraint(const std::vector<NSPgeometry::XYZ> &crd, const std::vector<double> &weights,double kres,
			double rmsdref, double rmsdswitch,int mode);
	StructRestraint(const std::vector<NSPgeometry::XYZ> &crd, const std::vector<std::pair<int,int>> &poslen,
			double kres,int mode=POSIMODE,double rmsdref=0.0);
	double posi_energy(const std::vector<NSPgeometry::XYZ> &crd,std::vector<NSPgeometry::XYZ> *forces) const;
	double rmsd2() const{return rmsd2_;}
	double totalrmsd_energy(const std::vector<NSPgeometry::XYZ> &crd, std::vector<NSPgeometry::XYZ> *forces) const;
	virtual double energy(const std::vector<NSPgeometry::XYZ> &crd,std::vector<NSPgeometry::XYZ> *forces) const{
		if(mode_==POSIMODE) ene_= posi_energy(crd,forces);
		else if(mode_==TOTALMODE) ene_= totalrmsd_energy(crd,forces);
		return ene_;
	}
	virtual void printinfo(std::ostream &os) const {
		os <<" RMSD restraint with kres, rmsdmin and rmsdswitch = "
				<< kres_/100 <<", " << rmsdref_*10 <<" and " <<rmsdswitch_*10<<std::endl;
	}
private:
    std::vector<NSPgeometry::XYZ> crdref_;
	std::vector<double> weights_;
	double wtot_{0.0};
	double kres_{0.0};
	double rmsdref_{0.0};
	double rmsdswitch_{10000.0};
	mutable double rmsd2_{0.0};
	int mode_{POSIMODE};
};
class RgRestraint:public ResTraint{
public:
	RgRestraint(double b=0,double k=0):rgbound_(b),kres_(k){;}
	RgRestraint(double b,double k,const std::vector<double> &w):rgbound_(b),kres_(k),w_(w){;}
	virtual double energy(const std::vector<NSPgeometry::XYZ> &crd, std::vector<NSPgeometry::XYZ> *forces) const {
// test
//		for(auto & f:*forces) f=NSPgeometry::XYZ(0,0,0);
		ene_= logrestr_energy(crd,forces);
/*/ test forces
		std::vector<NSPgeometry::XYZ> ctmp=crd;
		std::vector<NSPgeometry::XYZ> ftmp=*forces;
		for(int i=0;i<ctmp.size();++i){

			for(int m=0;m<3;++m){
				ctmp[i][m] += 0.0001;
	            double ep=logrestr_energy(ctmp,&ftmp);
	            ctmp[i][m] -=0.0002;
	            double em=logrestr_energy(ctmp,&ftmp);
	            std::cout << i<<" "<< m << (ep-em)/0.0002<<" "<<(*forces)[i][m] <<std::endl;
	            ctmp[i][m] +=0.0001;
			}
		}*/
		return ene_;
	}
	double logrestr_energy(const std::vector<NSPgeometry::XYZ> &crd, std::vector<NSPgeometry::XYZ> *forces) const;
/*	double energy1(const std::vector<NSPgeometry::XYZ> &crd, std::vector<NSPgeometry::XYZ> *forces
			,const std::vector<double> & w=std::vector<double>()) const;
	double energy2(const std::vector<NSPgeometry::XYZ> &crd, std::vector<NSPgeometry::XYZ> *forces
			,const std::vector<double> & w=std::vector<double>()) const;*/
	double &kres() {return kres_;}
	double kres() const {return kres_;}
	virtual void printinfo(std::ostream &os) const {
		os <<" Radius of gyration restraint with kres and rgbound = "
				<< kres_/2.4942 <<" and " <<rgbound_*10<<std::endl;
	}
private:
	double rgbound_;
	double kres_;
	std::vector<double>  w_;
};

struct DisRestraint:public ResTraint{
	int a1,a2;
	double kres;
	double r0;
	double r1;
	DisRestraint(int i=0,int j=0, double rij0=0,double rij1=0, double k=0):
		a1(i),a2(j),r0(rij0),r1(rij1),kres(k){;}
	virtual double energy(const std::vector<NSPgeometry::XYZ> &crd, std::vector<NSPgeometry::XYZ> *forces) const{
		if(r0>0) ene_= energy_attract(crd,forces);
		else ene_= energy_repulsion(crd,forces);
		return ene_;
	}
	double energy_attract(const std::vector<NSPgeometry::XYZ> &crd, std::vector<NSPgeometry::XYZ> *forces) const;
	double energy_repulsion(const std::vector<NSPgeometry::XYZ> &crd, std::vector<NSPgeometry::XYZ> *forces) const;
	virtual void printinfo(std::ostream &os) const {
		os <<" Distance restraint with atom1, atom2, kres, r0 and r1 = "
				<<a1<<", "<<a2<<", "<< kres <<", "<< r0 <<" and  "<< r1<<std::endl;
	}
};

/**
 * Distance restraint between the centers of two groups of atoms
 */
struct GrpDisRestraint:public ResTraint{
	std::vector<int> grp1,grp2;
	DisRestraint disres;
	GrpDisRestraint(){;}
	GrpDisRestraint(const std::vector<int> &g1, const std::vector<int> &g2,
			double rij0,double rij1,double k):
					grp1(g1),grp2(g2),disres(0,1,rij0,rij1,k){;}
	virtual double energy(const std::vector<NSPgeometry::XYZ> &crd, std::vector<NSPgeometry::XYZ> *forces) const;
	virtual void printinfo(std::ostream &os) const {
		os<<"Distance restraint between two groups of atoms with kres, r0 and r1 = "
				<<disres.kres <<", "<<disres.r0<<","<<disres.r1<<std::endl;
		os<<"Atoms in group 1: ";
		for(int i:grp1) os <<" "<<i;
		os<<std::endl;
		os<<"Atoms in group 2: ";
		for(int i:grp2) os <<" "<<i;
		os<<std::endl;
	}
};
struct HelixRestraint: public ResTraint{
//	std::vector<NSPintrct::SSRestraint> restraints;
	std::vector<DisRestraint> hbrestraints;
	const NSPintrct::IntrctMol *imol;
	virtual double energy(const std::vector<NSPgeometry::XYZ> &crd, std::vector<NSPgeometry::XYZ> *forces) const{
		ene_=0;
//		for(auto &res:restraints){
		for(auto &res:hbrestraints){
//			ene_+=res.energy(imol->sscodes,forces);
			ene_ +=res.energy(crd,forces);
		}
		return ene_;
	}
//	HelixRestraint(NSPintrct::IntrctMol *im, int flankingmcoff,double kres=500.0, double ptarget=0.95);
	HelixRestraint(NSPintrct::IntrctMol *im):imol(im){;}
	virtual void printinfo(std::ostream &os) const {
		os <<" Helix restraint with number of helix hbond distance restraints = "
				<<hbrestraints.size()<<std::endl;
	}
};
/*
 * restrain total number of contacts
 * the total number of  contact is calculated as the sum of the contact value or
 *  theta between two groups of atoms (see class GroupContact)
 */
struct ContactRestraint:ResTraint{
	std::vector<GroupContact> gcs;
	double nc0{0.0}; //target minimum number of contacts, if nc>nc0, restraint energy will be 0
	double kres{0.0};
	double resgdmin{ 0.5 };
	double resgdsmall{ 0.8 };
	double resgdoff{ 2.0 };
	ContactRestraint(){;}
	ContactRestraint(const NSPintrct::IntrctMol & imol,
			const NSPintrct::BlkSelector & iresidues_receptor,
			const NSPintrct::BlkSelector & iresidues_ligand,
			double nc0_i,double kres_i, double resgdmin_i, double resgdsmall_i, double resgdoff_i);
	virtual double energy(const std::vector<NSPgeometry::XYZ> &crd, std::vector<NSPgeometry::XYZ> *forces) const;
	virtual void printinfo(std::ostream &os) const{
		os <<"Contact Restraint with kres and ncontact0 =" << kres <<" and " <<nc0<<"; groupdistance with gdmin,gdsmall and gdoff ="<<resgdmin*10<<", "<<resgdsmall*10<<", "<<resgdoff*10<<std::endl;
		os<<"Atom group pairs for each potential contact: "<<std::endl;
		for(auto & g:gcs){
			 os<< "grp1:";
			 for( auto a:g.grpd.grps[0]) os <<" "<< a;
			 os <<"  grp2:";
			 for(auto a:g.grpd.grps[1]) os <<" "<<a;
			 os <<std::endl;
		}
	}
};
std::vector<std::shared_ptr<ResTraint> > readrestraints(NSPintrct::IntrctMol &imol,const std::string & resfile);
}

#endif /* RESTRAINTS_H_ */

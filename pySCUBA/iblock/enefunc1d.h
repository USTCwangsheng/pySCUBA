/*
 * enefunc1d.h
 *
 *  Created on: 2018年2月19日
 *      Author: hyliu
 */

#ifndef IBLOCK_ENEFUNC1D_H_
#define IBLOCK_ENEFUNC1D_H_
#include <cmath>
namespace NSPintrct{
/*
 * base class for one-dimensional math functions used in energy calculations
 */
class EneFunc1D{
public:
	/**
	 * x is the 1-D variable, dedx is delivered with the energy drivative
	 */
	virtual double energy(double x,double *dedr)=0;
	virtual double energyalt(double x, double *dedr){return 0.0;}
	virtual ~EneFunc1D(){;}
};

/**
 * no interaction
 */
class NullEne1D: public EneFunc1D{
public:
	virtual double energy(double x, double *dedr){
		*dedr=0.0;
		return 0.0;
	}
};

/**
 * a shifted Gaussian 1-D function
 */
class ShiftedG: public EneFunc1D {
public:
	ShiftedG(double r,double s, double e):x0_(r),sigma_(s),e0_(e){;}
	virtual double energy(double r, double *dedr){
		*dedr=0.0;
		if(r<=x0_) return e0_;
		double x=(r-x0_)/sigma_;
		if(x>4.5) return 0.0;
		double e=e0_*exp(-x*x);
		*dedr= -2.0*x*e/sigma_;
		return e;
	}
private:
	double x0_;
	double sigma_;
	double e0_;
};

/**
 * A 1-D function based on sine
 */
class SineEne: public EneFunc1D {
public:
	SineEne(double e0,double xp):e0_(e0),xperiod_(xp){;}
	virtual double energy(double r, double *dedr){
		*dedr=0.0;
		if(r>=xperiod_) return 0.0;
		double x=r/xperiod_*3.14159265;
		double e=e0_*sin(x);
		*dedr= e0_*cos(x)*3.14159265/xperiod_;
		return e;
	}
private:
	double e0_;
	double xperiod_;
};

/**
 * a 1-D function which is harmonic on the left side, and inverted gaussian on the right side
 *
 * This is used for sidechain atomic packing
 * e(r)= -emin_+[k_(r-rmin)]^m  if r<=rmin_
 *      -emin_*exp[-(r-rmin)^2/sigma_^2] if r>=rmin_
 *
 */
class HarmGauEne:public EneFunc1D{
public:
	HarmGauEne(double rmin,double emin,double lamda1,double lamda2,int m):
		rmin_(rmin),emin_(emin), m_(m)
	{	rs_=rmin_*lamda2;
		rcut_=rmin_+4.5*rs_;
		k_=1.0/pow((rmin_*lamda1),m);
		double tmp = pow(dedrconst_/(6.0*k_), 0.2);
		rconst_ = rmin_ - tmp;
		econst_ = k_ * pow(tmp, 6.0) - emin_;
	}
	virtual double energy(double r, double *dedr){
		*dedr=0.0;
		if(r>=rcut_) return 0.0;
		double e;
		if(r <= rconst_) {
		    double dr = r - rconst_;
		    e = econst_ - dr * dedrconst_;
		    *dedr = -1.0 * dedrconst_;
		} else if(r<=rmin_){
			double dr=r-rmin_;
			double drp=pow(dr,m_-1);
			e=-emin_+k_*drp*dr;
			*dedr=(double)(m_)*k_*drp;
			//if(*dedr<-4000.0) *dedr=-4000.0; //avoid too large repulsive forces
		} else {
			double x=(r-rmin_)/rs_;
			e=-emin_*exp(-x*x);
			*dedr=-e*2.0*x/rs_;
		}
		return e;
	}
	virtual double energyalt(double r, double *dedr){
		*dedr=0.0;
		if(r>rmin_) return 0.0;
		else return emin_+energy(r,dedr);
	}
private:
	double rmin_;
	double emin_;
	const double dedrconst_{500.0};
    double rconst_;
	double econst_;
	double k_;
	double rs_;
	int m_;
	double rcut_;
};

/**
 * 1-D energy function.
 *
 *  sigma_lj=rmin/(2^1/6)
 * sigma_g=rmin_*lamda2
 * e(r)= 4e0(r12/sigma_lj^12-r6/sigma_lj^6)+e0-emin_  if r<=rmin)
 *      -emin_*exp[-(r-rmin)^2/sigma_g^2] if r>=rmin_
*/
class LJGauEne:public EneFunc1D{
public:
	LJGauEne(double rmin,double emin,double lamda2):
		rmin_(rmin),emin_(emin)
	{	sigma_g_=rmin_*lamda2;
		sigma_lj_=rmin_/pow(2.0,1/6.0);
		sigma_lj6_=pow(sigma_lj_,6);
		sigma_lj12_=sigma_lj6_*sigma_lj6_;
		rcut_=rmin_+4.5*sigma_g_;
		rconst_=sigma_lj_*0.85;
		double r6=pow(rconst_,6);
		econst_=4.0*e0_*(sigma_lj12_/(r6*r6)-sigma_lj6_/r6)+e0_-emin_;
		dedrconst_=4.0*e0_*(-12.0*sigma_lj12_/(r6*r6*rconst_)
			+6.0*sigma_lj6_/(r6*rconst_));
	}
	virtual double energy(double r, double *dedr){
		*dedr=0.0;
		if(r>=rcut_) return 0.0;
		double e;
		if(r <= rconst_) {
		    double dr = r - rconst_;
		    e = econst_ + dr * dedrconst_;
		    *dedr = dedrconst_;
		} else if(r<rmin_){
			double r3=r*r*r;
			double r6=r3*r3;
			e=4.0*e0_*(sigma_lj12_/(r6*r6)-sigma_lj6_/r6)+e0_-emin_;
			*dedr=4.0*e0_*(-12.0*sigma_lj12_/(r6*r6*r)
				+6.0*sigma_lj6_/(r6*r));
			//if(*dedr<-4000.0) *dedr=-4000.0; //avoid too large repulsive forces
		} else {
			double x=(r-rmin_)/sigma_g_;
			e=-emin_*exp(-x*x);
			*dedr=-e*2.0*x/sigma_g_;
		}
		return e;
	}
	virtual double energyalt(double r, double *dedr){
		*dedr=0.0;
		if(r>rmin_) return 0.0;
		else return emin_+energy(r,dedr);
	}
private:
	double rcut_;
	double rmin_;
	double sigma_lj_;
	double sigma_lj12_;
	double sigma_lj6_;
	double sigma_g_;
	double emin_;
	double e0_{0.4};
	double rconst_;
	double econst_;
	double dedrconst_;
};

/**
 * 1-D energy function
 */
class SwitchGauEne:public EneFunc1D{
public:
	SwitchGauEne(double rmin,double emin,double rcore,double emax,double lamda2):
		rmin_(rmin),emin_(emin), rcore_(rcore),emax_(emax)
	{	dr_=rmin_-rcore_;
		de_=emax_-emin_;
		rs_=rmin_*lamda2;
		rcut_=rmin_+4.5*rs_;
	}
	virtual double energy(double r, double *dedr){
		*dedr=0.0;
		if(r>=rcut_) return 0.0;
		double e;
		if(r<=rmin_){
			if(r <= rcore_) {
				*dedr = 0.0;
				e = emax_;
			} else {
				double x = (r - rcore_) / dr_;
				if (x <= 0.5) {
					e = emax_-2.0 * x * x*de_;
					*dedr = -4.0 * x*de_ / dr_;
				} else {
					e = emin_+ 2.0 * (1 - x) * (1 - x)*de_;
					*dedr = -4.0 * (1 - x)*de_ / dr_;
				}
			}
		} else {
			double x=(r-rmin_)/rs_;
			e=emin_*exp(-x*x);
			*dedr=-e*2.0*x/rs_;
		}
		return e;
	}
private:
	double rmin_;
	double emin_;
	double emax_;
	double rcore_;
	double rs_;
	double de_;
	double dr_;
	double rcut_;
};
}



#endif /* IBLOCK_ENEFUNC1D_H_ */

/*
 * bondedenergies.cpp
 *
 *  Created on: 2019年12月2日
 *      Author: hyiu
 */
#include "iblock/intrctblck.h"
#include "geometry/calculators.h"
#include "iblock/backboneff.h"
#include "iblock/vsctypes.h"

using namespace NSPintrct;
template<>
double Bond::energy(const std::array<NSPgeometry::XYZ, 2> &crd,
		std::array<NSPgeometry::XYZ, 2>* dedx) const {
	std::vector<NSPgeometry::XYZ> deriv;
	double b = NSPgeometry::distance(crd.at(0), crd.at(1), &deriv);
	double db = b - v0;
	double e = 0.5 * k * db * db;
	double kbdb = k * db;
	(*dedx)[0] = kbdb * deriv[0];
	(*dedx)[1] = kbdb * deriv[1];
	return e;
}
template<>
double Angle::energy(const std::array<NSPgeometry::XYZ, 3> &crd,
		std::array<NSPgeometry::XYZ, 3>* dedx) const {
	std::vector<NSPgeometry::XYZ> deriv;
	double ct = cos_angle(crd[0], crd[1], crd[2], &deriv);
	double dct = ct - v0;
	double e = 0.5 * k * dct * dct;
	double ktdct = k * dct;
	for (int i = 0; i < 3; ++i)
		(*dedx)[i] = ktdct * deriv[i];
	return e;
}
template<>
double ImpDih::energy(const std::array<NSPgeometry::XYZ, 4> &crd,
		std::array<NSPgeometry::XYZ, 4>* dedx) const {
	std::vector<NSPgeometry::XYZ> deriv;
	double p = torsion(crd[0], crd[1], crd[2], crd[3], &deriv);
	double dp = p - v0;
	if (dp > 3.14159265)
		dp -= 2 * 3.14159265;
	else if (dp < -3.14159265)
		dp += 2 * 3.14159265;
	double e = 0.0;
	double kf; // = k * dp;
	if (dp > 0.09)
	{
		e = 0.5 * k * 0.0081 + k * 0.09 * (dp-0.09);
		kf = k * 0.09;
	}
	else if(dp < -0.09)
	{
		e = 0.5 * k * 0.0081 - k * 0.09 * (dp+0.09);
		kf = -k * 0.09;
	}
	else
	{
		e = 0.5 * k * dp * dp;
		kf = k * dp;
	}

/*	for (auto c:crd){
		std::cout <<c.toString()<<" "<<p*180/3.14159265<<std::endl;
		std::cout<<angle(crd[0],crd[1],crd[2]) <<"    "<<angle(crd[1],crd[2],crd[3]) <<std::endl;
	}
	std::cout <<"IMP:"<< dp*180/3.14159265<<std::endl;*/
	for (int i = 0; i < 4; ++i)
		(*dedx)[i] = kf * deriv[i];
	return e;
}
double AtomIP::pairenergy(const AtomIP &sa1, const NSPgeometry::XYZ &pos1,
		const AtomIP&sa2, const NSPgeometry::XYZ &pos2,
		std::array<NSPgeometry::XYZ, 2>*dedx) {
	static const double fac = pow(2, 1 / 6.0);
	static const double DEDRMAX{400.0};
	static const double rcut2_m = 0.25; // in nm**2
	NSPgeometry::XYZ dr = pos2 - pos1;
	double dr2=dr.squarednorm();
	 if ( dr2> rcut2_m) {
			return 0.0;
	 }
	dedx->fill(NSPgeometry::XYZ { 0.0, 0.0, 0.0 });
	double sigma;
	if ((sa1.hbdonor && sa2.hbacceptor) || (sa1.hbacceptor && sa2.hbdonor)) {
		sigma = 0.5 * (sa1.sigmahb + sa2.sigmahb);
	} else if ((sa1.hbdonor && sa2.hbdonor)
			|| (sa1.hbacceptor && sa2.hbacceptor)) {
		sigma = 0.55 * (sa1.sigma + sa2.sigma);
	} else {
		sigma = 0.5 * (sa1.sigma + sa2.sigma);
	}
	double rcut2 = sigma * sigma * 1.2;
	if (dr2> rcut2) {
		return 0.0;
	}
	double eps = sqrt(sa1.eps * sa2.eps);
	double sl = fac * sigma;
	if (dr.squarednorm() > sl * sl) {
		return 0.0;
	}
	std::vector<NSPgeometry::XYZ> deriv;
	double r = NSPgeometry::distance(pos1, pos2, &deriv);
	double r2 = r * r;
	double r3 = r2 * r;
	double r5 = r2 * r3;
	double r6 = r3 * r3;
	double r7 = r * r6;
	double r12 = r6 * r6;
	double r13 = r * r12;
	double sigma3 = sigma * sigma * sigma;
	double sigma6 = sigma3 * sigma3;
	double sigma12 = sigma6 * sigma6;
	double ene = 4.0 * eps * (sigma12 / r12 - sigma6 / r6) + eps;
	double dedr = 4.0 * eps * (-12.0 * sigma12 / r13 + 6.0 * sigma6 / r7);
	if (dedr < - DEDRMAX)
		dedr = -DEDRMAX;
	(*dedx)[0] = dedr * deriv[0];
	(*dedx)[1] = dedr * deriv[1];
	return ene;
}
PrePro::PrePro(const BlckTopo &topo){
	static BackBoneFF backboneff;
	std::array<AtomIdx3D,3> aidx;
	aidx[0]=AtomIdx3D{0,0,topo.natoms-2};
	aidx[1]=AtomIdx3D{0,1,0};
	aidx[2]=AtomIdx3D{0,1,BlckTopo::CD_PRO};
	angles.push_back(Angle(aidx,cos(backboneff.t0_ocn*3.14159265/180.0),
			2 * KBT * KANG_FAC * backboneff.kt_ocn));
	std::array<AtomIdx3D,4> pidx;
	pidx[0]=AtomIdx3D{0,0,1};
	pidx[1]=AtomIdx3D{0,0,topo.natoms-2};
	pidx[2]=AtomIdx3D{0,1,0};
	pidx[3]=AtomIdx3D{0,1,BlckTopo::CD_PRO};
	impdihs.push_back(ImpDih(pidx,0.0,
			2 * KBT * KANG_FAC * backboneff.kp_cacnca));
	pidx[0]=AtomIdx3D{0,0,topo.natoms-1};
	pidx[1]=AtomIdx3D{0,0,topo.natoms-2};
	pidx[2]=AtomIdx3D{0,1,0};
	pidx[3]=AtomIdx3D{0,1,BlckTopo::CD_PRO};
	impdihs.push_back(ImpDih(pidx,3.14159265,
			2 * KBT * KANG_FAC * backboneff.kp_cacnca));
}
void PrePro::extraforces(const IntrctBlck &blck,const IntrctPara &param,
		std::vector<NSPgeometry::XYZ> *dedx_sys){
	PrePro prepro(*blck.topo);
	blck.energies[IntrctBlck::ANGLE] += blck.bondedforce(param,
			prepro.angles,
			dedx_sys);
	if(blck.precis){
		prepro.impdihs[0].v0=3.14159265;
		prepro.impdihs[1].v0=0.0;
	}
	blck.energies[IntrctBlck::IMPDIH] += blck.bondedforce(param,
			prepro.impdihs,
			dedx_sys);
}
double IntrctBlck::blckstericforces(const IntrctPara &param,
		const IntrctBlck &blk2, const AtomNbList &nblist,
		std::vector<NSPgeometry::XYZ> *dedx_sys) const {
	double esum = 0.0;
	assert(nblist.size() == natoms());
	double w = param.weight_steric;
	std::array<NSPgeometry::XYZ, 2> dedx;
	for (int i = 0; i < natoms(); ++i) {
		int iidx = aoffset + i;
		for (auto j : nblist.at(i)) {
			double ene =w*AtomIP::pairenergy(topo->atomips[i], crds[i],
							blk2.topo->atomips[j], blk2.getcrd(j), &dedx);
			esum+=ene;
			int jidx = blk2.aoffset + j;
           if(param.enedetails){
        	   Results::ofstreams(param.jobname,"MCStericEne")
        			<<atomstring(i)<<"-"<<blk2.atomstring(j)<<" e = "<<ene<<std::endl;
           }
			(*dedx_sys)[iidx] = (*dedx_sys)[iidx] + w * dedx[0];
			(*dedx_sys)[jidx] = (*dedx_sys)[jidx] + w * dedx[1];
		}
	}
	return esum;
}

double IntrctBlck::stericforces(const IntrctPara &param,
		std::vector<NSPgeometry::XYZ> *dedx_sys) const {
	double esum = 0.0;
	for (auto & entry : nblist_M) {
		double ene= blckstericforces(param, (*molsystm_)[entry.first], entry.second,
				dedx_sys);
		this->energies[STERIC]+=0.5*ene;
		(*molsystm_)[entry.first].energies[STERIC]+=0.5*ene;
		esum+=ene;
	}
	for(auto &entry:nblist_L){
		double ene= blckstericforces(param, (*molsystm_)[entry.first], entry.second,
				dedx_sys);
		this->energies[STERIC]+=0.5*ene;
		(*molsystm_)[entry.first].energies[STERIC]+=0.5*ene;
		esum+=ene;
	}
	return esum;
}
#include "iblock/intrctmol.h"
void IntrctMol::forces_cov(const IntrctPara &param,std::vector<NSPgeometry::XYZ> &dedx) const{
	assert(offsetsok_);
	auto &chains=mol_.D_;
	for(auto &c:chains){
		for(auto &r:c) {
			r.energies[IntrctBlck::BOND]=0.0;
			r.energies[IntrctBlck::ANGLE]=0.0;
			r.energies[IntrctBlck::IMPDIH]=0.0;
		}
	}
	for(auto &c:chains){
		for(auto &r:c){
//			std::cout <<"residue: "<<residue++<<std::endl;
			if(r.activemod==IntrctBlck::INACTIVE &&
					!mainchainwindowactive(r.getidx2d(),1)) continue;
			r.covalentforces(param,&dedx);
		}
	}
}
void IntrctMol::forces_steric(const IntrctPara &param,std::vector<NSPgeometry::XYZ> &dedx) const{
	assert(offsetsok_);
	auto &chains=mol_.D_;
	int residue=0;
	for(auto &c:chains){
		for(auto &r:c) r.energies[IntrctBlck::STERIC]=0.0;
	}
	for(auto &c:chains){
		for(auto &r:c){
//			std::cout <<"residue: "<<residue++<<std::endl;
			r.stericforces(param, &dedx);
		}
	}
}

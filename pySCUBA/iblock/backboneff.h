/*
 * backboneff.h
 *
 *  Created on: 2018年1月25日
 *      Author: hyliu
 */

#ifndef IBLOCK_BACKBONEFF_H_
#define IBLOCK_BACKBONEFF_H_
#define NATOMS_PER_BACKBONESITE 4
namespace NSPintrct{
/**
 * forcefield-like interaction parameters for peptide backbone
 *
 * values like bond lengths, bond angles, L-J parameters, force constants are provided here
 *
 * Energy unit is KBT.
 */
struct BackBoneFF {
	int natomspersites { NATOMS_PER_BACKBONESITE };
	double b0_nca { 1.460 };
	double kb_nca { 1000.0 };
	double b0_cac { 1.526 };
	double kb_cac { 1000.0 };
	double b0_co { 1.233 };
	double kb_co { 1000.0 };
	double b0_cn { 1.330 };
	double kb_cn { 1000.0 };
	double t0_ncac { 111.04 };
	double kt_ncac { 0.01417 };
	double t0_caco { 120.54 };
	double kt_caco { 0.1232 };
	double t0_cacn { 116.53 };
	double kt_cacn { 0.05089};
	double t0_ocn { 122.90 };
	double kt_ocn { 0.126 };
	double t0_cnca { 121.52 };
	double kt_cnca { 0.08};
	double p_cncao { 0.0 };
	double kp_cncao { 0.06 };
	double p_cacnca { 180.0 };
	double kp_cacnca { 0.05};
	double p_ocnca { 0.0 };
	double kp_ocnca { 0.01 };
	double sigma_n { 3.0 };
	double sigma_nhb { 2.8 };
	double sigma_ca { 3.1 };
	double sigma_c { 3.0 };
	double sigma_o { 3.0 };
	double sigma_ohb { 2.8 };
	double eps { 0.40 };
};

}
#endif /* SD_BACKBONEFF_H_ */

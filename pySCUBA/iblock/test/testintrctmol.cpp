/*
 * testintrctmol.cpp
 *
 *  Created on: 2019年12月4日
 *      Author: hyiu
 */
#include "iblock/intrctmol.h"
#include "iblock/molsystmpara.h"
//#include "fullsite/fullsite.h"
#include <string>
using namespace NSPproteinrep;
using namespace NSPintrct;
/**
 * This program illustrate  how to use the IntrctMol class to
 *  calculate SCUBA energies and atomic forces of a protein
 *
 * Usage:
 * testintrctmol controling_parameter_filename [jobname]
 * The first command line argument controling_parameter_filename gives the file
 * from which most input parameters are read.
 * The second command line argument jobname is optional. If given, the output files
 * will be named using this given jobname. Otherwise, they will be named using the
 * "jobname" defined in the controlling parameter file.
 *
 */
int main(int argc, char** argv){
	std::string parfile(argv[1]);
	std::string jobname("auto");
	if(argc>2) jobname=std::string(argv[2]); //jobname defined by command line argument
	NSPdataio::ControlFile cf;
	cf.readfile(parfile);  //read controlling paramters into memory
	MolSystmPara molpara=makemolsystmparam(std::string(),cf);  //extract parameters on how to define the molecule system
	if(jobname=="auto"){
		jobname=molpara.jobname;  //jobname from controlling parameter
	} else {
		molpara.jobname=jobname; //replace jobname in molecular system parameter from command line
	}
	IntrctPara para_o=makeintrctparam(jobname,cf); // extract parameters to define the energy function
	std::shared_ptr<IntrctMol> imol=make_molsystm(molpara,para_o);  //make the molecule system data structure
	std::vector<NSPgeometry::XYZ> dedx;  //variable for storing the derivatives of energy with respect to coordinates
	std::array<double,IntrctBlck::ENESIZE> energies; //variable for storing energy components
	imol->forces_all(dedx);  //calculate the energy (components stored in imol ) and derivatives (delivered into dedx)
/* uncomment the following part to compare the analytically calculate energy derivatives with
// finite difference derivatives
    std::array<NSPgeometry::XYZ,3>  dr{NSPgeometry::XYZ(0.00001,0.0,0.0),
		NSPgeometry::XYZ(0.0,0.00001,0.0),NSPgeometry::XYZ(0.0,0.0,0.00001)}; // coordinate changes to finite difference evaluation of derivatives
    para_o.enedetails=false;
	para_o.calc_sscodes=false;
	para_o.calc_nblists=false;
	std::vector<NSPgeometry::XYZ>dedxtmp;
	std::vector<int> atoms{623, 624, 631,632,633,634,638,639};
	for(int i:atoms){
//	for(int i=0;i<mol.natoms();++i){
		std::cout <<"atom "<<i<<" "<<mol.atomname(i)<<std::endl;
		for(int m=0;m<3;++m){
			crds[i]=crds[i]+dr[m];
			mol.forces_all(para_o,crds,dedxtmp);
			double e1=mol.sum_energies(energies);
			crds[i]=crds[i]-2.0*dr[m];
			mol.forces_all(para_o,crds,dedxtmp);
			double e2=mol.sum_energies(energies);
			std::cout << dedx[i][m]<<" "<< (e1-e2)/0.00002<<std::endl;
			crds[i]=crds[i]+dr[m];
		}
	}*/
}



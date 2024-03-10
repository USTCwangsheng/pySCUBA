/*
 * testaaconformer.cpp
 *
 *  Created on: 2020年7月17日
 *      Author: hyiu
 */
#include "proteinrep/aaconformer.h"
#include "proteinrep/rotamerutil.h"
#include "geometry/rigidtransform.h"
using namespace NSPproteinrep;
using namespace NSPdesignseq;
using namespace NSPgeometry;

int main(){
	AAConformer his=make_aaconformer(randomrotamer("HIS"));
	AAConformer asp=make_aaconformer(randomrotamer("ASP"));
	AAConformer ser=make_aaconformer(randomrotamer("SER"));
	std::vector<XYZ> joints1;
	joints1.push_back(his.getglobalcrd()["NE2"]);
	joints1.push_back(his.getglobalcrd()["CE1"]);
	joints1.push_back(his.getglobalcrd()["ND1"]);
	joints1.push_back(ser.getglobalcrd()["OG"]);
	joints1.push_back(his.getglobalcrd()["CB"]);
	joints1.push_back(his.getglobalcrd()["CA"]);
	std::vector<double> intcrds1{3.14159265,2.09,3.0,3.14159265,2.1,3.14159265};
	RigidTransform rt1(joints1,intcrds1);
	for(auto &x:ser.getglobalcrd()) rt1.apply(&x.second);

	std::vector<std::string> aspanchors{"OD1","CG","CB"};
	std::map<std::string,XYZ> asprefcrds;
	asprefcrds["OD1"]=InternaltoXYZ(his.getglobalcrd()["NE2"], his.getglobalcrd()["CE1"],
			his.getglobalcrd()["ND1"],2.8, 120.0*3.14159265/180.0,3.14159265);
	asprefcrds["CG"]=InternaltoXYZ(asprefcrds["OD1"],
			his.getglobalcrd()["NE2"], his.getglobalcrd()["CE1"],
			1.25, 150.0*3.14159265/180.0,3.14159265);
	asprefcrds["CB"]=InternaltoXYZ(asprefcrds["CG"],asprefcrds["OD1"],
				his.getglobalcrd()["NE2"],
				1.50, 119.0*3.14159265/180.0,0.5*3.14159265);
	asp.moveglobalcrd(asprefcrds,aspanchors);
	int natoms=1;
	auto hispdbrcds=his.make_pdbrecords('A',1,' ',natoms);
	natoms+=hispdbrcds.size();
	auto serpdbrcds=ser.make_pdbrecords('A',2,' ',natoms);
	natoms +=serpdbrcds.size();
	auto asppdbrcds=asp.make_pdbrecords('A',3,' ',natoms);
    for(auto &r:hispdbrcds) std::cout<<r.toString()<<std::endl;
    for(auto &r:serpdbrcds) std::cout<<r.toString()<<std::endl;
    for(auto &r:asppdbrcds) std::cout<<r.toString()<<std::endl;
}




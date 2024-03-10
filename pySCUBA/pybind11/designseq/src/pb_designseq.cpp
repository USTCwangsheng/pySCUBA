/*
 * @Author: ws
 * @2022-05-23
 */

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>

#include "designseq/ProteinRep.h"
#include "designseq/StructureInfo.h"
#include "designseq/S1EnergyTable.h"
#include "geometry/xyz.h"
#include "designseq/AtomLib.h"
#include "designseq/StructureInfo.h"
#include "designseq/S2MatrixFinder.h"
#include "designseq/AtomicEnergyCalcular.h"
#include "designseq/DesignTemplate.h"
#include "designseq/AAProbabilityArray.h"
#include "designseq/SeqMinimize.h"
#include "designseq/CmdArgs.h"

#include "designseq/test/TestDesignSeq.cpp"
#include "dstl/randomengine.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <time.h>

using namespace std;
namespace py = pybind11;
using namespace NSPdesignseq;
using namespace NSPgeometry;
namespace desi = NSPdesignseq;
namespace geo = NSPgeometry;

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;



void SCUBAdesignseq(std::string& pdbname,std::string& outname, std::string& logfile, int & outn, float& vdw,
	double & div, double& nat,std::string& inresfile ,  std::string& inparafile, std::string& aapropfile) {
	/*
 * input is PDB structure
 */

	int n = 1;  //n个结构
	n = outn;

	clock_t start = clock();   //开始计时


	string s = pdbname;//输入的pdb
	string pdbID = "unk";
	string output = outname;  //输出名
	string logFile = logfile;
	if (output.length() < 4 || output.substr(output.length() - 4, 4) != ".pdb")  //处理输出文件名后缀pdb
		output = output + ".pdb";


	PDB pdb(s, pdbID);
	cout << "start" << endl;

	std::string ssseq;
	vector<BackBoneSite*> bsList;
	pdb2BackboneSiteList(&pdb, bsList, ssseq);
	// chainid of bsLiist is the same with number in pdbfile
	// resseq is sequence reranked

	std::vector<AANumberControlInMC> aacmc;
	string aapropFile = aapropfile;
	if (aapropFile != "") {
		readaaprop(aapropFile, ssseq, bsList, aacmc);
	}

	DesignParameters* dp;
	S1EnergyTable s1Etable;
	S2MatrixFinder s2Etable;

	DesignTemplate* dt;

	string parafile = "";
	parafile = inparafile;
	if (parafile != "") {
		dp = new DesignParameters(parafile);
	}
	else {
		dp = new DesignParameters();
	}
	cout << "start initializing design template" << endl;

	double vdwres = 1.0;
	vdwres = vdw;
	cout << "vdwres: " << vdwres << endl;

	string resfile = "";
	resfile = inresfile;

	if (resfile != "") {
		dt = new DesignTemplate(bsList, resfile, dp, s1Etable, s2Etable, vdwres);
	}
	else {
		dt = new DesignTemplate(bsList,dp, s1Etable, s2Etable, vdwres);
	}

	string natSeq = "";
	ResName rn;
	for (int i = 0;i < bsList.size();i++) {
		char c = rn.triToSin(bsList[i]->resname);
		natSeq = natSeq + c;
	}
	cout << "start loading S1S2" << endl;

	dt->loadS1S2(s1Etable, s2Etable);
	cout << "start loading single Residue Energy" << endl;
	dt->loadSingleResidueEnergy();
	cout << "start loading pairwise energy" << endl;
	dt->loadPairwiseEnergy();
	cout << "start MC" << endl;

	DesignMC mc(dt);
	RotSequence unit(dt->resNum);
	ofstream out(logFile, ios::out);

	double divCutoff = 1.0;   //默认应是1
	double sigmaDiv = 0.02;

	double eNatSeq = 0;

	double natIdt = -100;    //默认应是-100不执行
	double natSeqWt = 0;

	divCutoff = div;
	if(divCutoff!=1){ 
		cout << "divCutoff: " << divCutoff << endl; 
	}
	natIdt = nat;


	if (natIdt != -100) {
		if (natIdt < 0.0 || natIdt > 1.0) {
			cout << "-nat : specify the sequence identity between native template and your design result, this value should be between 0 and 1 " << endl;
			exit(1);
		}
		natSeqWt = 5.0 * (0.37 - natIdt);
		out << "native sequence identity energy: " << natSeqWt << endl;
		dt->updateNativeSequenceEnergy(natSeqWt);
	}

	vector<string> preSeqs;
	string preSeq = "";
	for (int i = 1;i <= n;i++)
	{
		char x[10];
		char s[1000];
		sprintf(x, "000%d", i);
		string xx(x);
		xx = xx.substr(xx.length() - 3, 3);
		string fileName = output;
		if (n > 1)
			fileName = fileName.substr(0, output.length() - 4) + "-" + xx + ".pdb";
		double score;
		if (preSeqs.size() > 0 && natIdt!=-100) {
			double preIdt = sequenceIdentity(natSeq, preSeq);
			if (abs(preIdt - natIdt) > 0.02) {
				natSeqWt = 5.0 * (preIdt - natIdt);
				dt->updateNativeSequenceEnergy(natSeqWt);
			}
		}

		if (divCutoff!=1) {
			if (aapropFile!="")
				score = mc.mcRunWithIdtAACountRestrain(&unit, preSeqs, divCutoff, sigmaDiv, aacmc);
			else
				score = mc.mcRunWithIdtRestrain(&unit, preSeqs, divCutoff, sigmaDiv);
		}
		else {
			if (aapropFile != "")
				score = mc.mcRunWithAACountRestrain(&unit, aacmc);
			else
				score = mc.mcRun(&unit);
		}



		double Es1, Erot, Eref;
		unit.calcEnergyComponents(dt, &Es1, &Erot, &Eref);

		string seq = unit.toAASequence();
		preSeqs.push_back(seq);
		preSeq = seq;
		sprintf(s, "%s %10.3f %10.3f", seq.c_str(), score, score - Eref);
		string result(s);
		out << result << endl;
		mc.printPDB(&unit, fileName);
	}

	out.close();

	delete dt;

	clock_t end = clock();
	cout << "time: " << (float)(end - start) / CLOCKS_PER_SEC << endl;

}


PYBIND11_MODULE(pb_designseq, m) {
	m.def("readaaprop", &readaaprop);
	m.def("SCUBAdesignseq", &SCUBAdesignseq);


};
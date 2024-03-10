/*
 * TestDesignSeq.cpp
 *
 *  Created on: 2017��12��27��
 *      Author: notxp
 */


#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <time.h>
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
#include <sstream>

using namespace std;
using namespace NSPdesignseq;
using namespace NSPgeometry;

void readaaprop(std::string parfile, std::string ssseq, vector<BackBoneSite*> &bsList,
        std::vector<AANumberControlInMC> &aacmc) {
    std::map<std::string,char> aa31 {
        {"GLY",'G'},
        {"ALA",'A'},
        {"VAL",'V'},
        {"LEU",'L'},
        {"ILE",'I'},
        {"SER",'S'},
        {"THR",'T'},
        {"CYS",'C'},
        {"MET",'M'},
        {"ASP",'D'},
        {"GLU",'E'},
        {"ASN",'N'},
        {"GLN",'Q'},
        {"LYS",'K'},
        {"ARG",'R'},
        {"HIS",'H'},
        {"PRO",'P'},
        {"PHE",'F'},
        {"TYR",'Y'},
        {"TRP",'W'}
    };

    std::map<char,std::set<int>> ssn{{'H',std::set<int>()},{'E',std::set<int>()},{'C',std::set<int>()}};
    for(int i=0;i<ssseq.size();i++) {
        char c=ssseq[i];
        ssn.at(c).insert(i);
    }
    std::cout <<"Number of Site in Helix: " <<ssn.at('H').size() <<std::endl;
    std::cout <<"Number of Site in Strand: " <<ssn.at('E').size() <<std::endl;
    std::cout <<"Number of Site in Loop: " <<ssn.at('C').size() <<std::endl;
    std::map<std::pair<char,int>,int> ids;
    for(int i=0;i<bsList.size();i++) {
        ids.insert({{bsList[i]->chainid,bsList[i]->resseq},i});
    
    }

    int minnum=5;
    int maxaanum=2;

    std::string readline;
    std::stringstream ss;
    std::ifstream ifs(parfile);
    while(std::getline(ifs,readline)) {
        if(readline.empty()) continue;
        if(readline[0]=='#') continue;
        std::cout <<readline <<std::endl;

        ss << readline;
        std::string a3;
        ss >> a3;
        char a1=aa31.at(a3);

        int nr;
        ss >> nr;

        std::vector<std::pair<char,int>> rid;
        for(int i=0;i<nr;i++) {
            char cid;
            int st,en;
            ss >> cid >> st >> en;
            for(int i=st;i<=en;i++) {
                std::pair<char,int> p{cid,i};
                if(ids.find(p)==ids.end()) {
                    std::cout <<"\tAAProp: Can Not Find Residue: " <<cid <<' ' <<i <<std::endl;
                    continue;
                }
                rid.push_back(p);
            }
        }

        char sslb;
        ss >> sslb;

        double min, max;
        ss >> min >> max;

        ss.clear();

        std::set<int> si;
        for(auto &p:rid) {
            int ix=ids.at(p);
            if(sslb=='A') si.insert(ix);
            else {
                char s=ssseq[ix];
                if(s==sslb) si.insert(ix);
            }
        }

        if(si.size()<minnum) continue;
        AANumberControlInMC ac;
        ac.sites = si;
        ac.aa = a1;
        ac.max = (int)(max*(double)(si.size()));
        if(ac.max<maxaanum) ac.max = maxaanum;
        ac.min = (int)(min*(double)(si.size()));
        //if (min > 0.1 && ac.min == 0) ac.min = 1;
        ac.num = 0;
        aacmc.push_back(ac);
    }
    ifs.close();

    std::cout <<"Amino Acid Number Restrain:" <<std::endl;
    for(auto &ac:aacmc) {
        std::cout <<ac.aa <<'\t' <<ac.min <<'\t' <<ac.max <<'\t' <<ac.sites.size() <<std::endl;
        for(auto &i:ac.sites) std::cout <<'\t' <<i;
        std::cout <<std::endl;
    }
}

double sequenceIdentity(const string& seqA, const string& seqB){
    if(seqA.length() != seqB.length())
        return 0.0;
    int sam = 0;
    for(int i=0;i<seqA.length();i++){
        if(seqA[i] == seqB[i])
            sam ++;
    }
    return 1.0*sam/seqA.length();
}

int main(int argc, char** args){
    /*
     * input is PDB structure
     */
    CmdArgs cmd(argc, args);

    if(!cmd.specifiedOption("-in") || !cmd.specifiedOption("-out") || !cmd.specifiedOption("-log")){
        cout << "DesignSeq -in $INPUTPDB -out $OUTPUTPDB -log $LOGFILE" << endl;
        cout << endl;
        cout << "optional tags:" << endl;
        cout << "-n: specify the sequence number" << endl;
        cout << "-div: specify the sequence diversity between designed results" << endl;
        cout << "-nat: specify the sequence identity between native sequence and designed sequence" << endl;
        cout << "-resFile: specify the amino acid choice on each position" << endl;
        cout << "-para: specify the parameter file" << endl;
        cout << "-aaprop: specify number of times for amino acid in the specified region" << endl;
        cout << "-VDW: Rescale r1+r2 when less than r0" << endl;
        cout << endl;
        return 1;
    }



    int n = 1;
    if(cmd.specifiedOption("-n")){
        n = atoi(cmd.getValue("-n").c_str());
    }


    clock_t start = clock();


    string s = cmd.getValue("-in");
    string pdbID = "unk";
    string output = cmd.getValue("-out");
    string logFile = cmd.getValue("-log");

    if(output.length() < 4 || output.substr(output.length()-4, 4) != ".pdb")
        output = output + ".pdb";



    PDB pdb(s, pdbID);

    cout << "start" << endl;

//    ProteinChain* pc = pdb.getFirstChain();
//    string natSeq = pc->getSequence();


//    proteinchain2BackboneSiteList(pc, bsList0);

    std::string ssseq;
    vector<BackBoneSite*> bsList;
    pdb2BackboneSiteList(&pdb, bsList, ssseq);
    // chainid of bsLiist is the same with number in pdbfile
    // resseq is sequence reranked

    std::vector<AANumberControlInMC> aacmc;
    if(cmd.specifiedOption("-aaprop")) {
        string aapropFile = cmd.getValue("-aaprop");
        readaaprop(aapropFile,ssseq,bsList,aacmc);
    }

    DesignParameters* dp;
    S1EnergyTable s1Etable;
    S2MatrixFinder s2Etable;

    DesignTemplate* dt;

    if(cmd.specifiedOption("-para")) {
        string paraFile = cmd.getValue("-para");
        dp = new DesignParameters(paraFile);
    }
    else
        dp = new DesignParameters();

    cout << "start initializing design template" << endl;

    double vdwres = 1.0;
    if (cmd.specifiedOption("-VDW")) {
        string srvdw = cmd.getValue("-VDW");
        float rvdw = std::stod(srvdw);
        vdwres = rvdw;
        cout << "vdwres: "<< vdwres << endl;
    }
    if(cmd.specifiedOption("-resFile")){
        string resFile = cmd.getValue("-resFile");
        dt = new DesignTemplate(bsList, resFile, dp, s1Etable, s2Etable, vdwres);
    }
    else{
        dt = new DesignTemplate(bsList, dp, s1Etable, s2Etable, vdwres);
    }



    string natSeq = "";
    ResName rn;
    for(int i=0;i<bsList.size();i++){
        char c = rn.triToSin(bsList[i]->resname);
        natSeq = natSeq + c;
    }
    //    cout << "vdwres: " << vdwres << endl;
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

    double divCutoff = 1.0;
    double sigmaDiv = 0.02;

    double eNatSeq = 0;

    double natIdt = -1;
    double natSeqWt = 0;

    if(cmd.specifiedOption("-div"))
        divCutoff = atof(cmd.getValue("-div").c_str());
    cout << "divCutoff: " << divCutoff << endl;
    if(cmd.specifiedOption("-nat"))
    {
        natIdt = atof(cmd.getValue("-nat").c_str());
        if(natIdt < 0.0 || natIdt > 1.0) {
            cout << "-nat : specify the sequence identity between native template and your design result, this value should be between 0 and 1 " << endl;
            exit(1);
        }
        natSeqWt = 5.0*(0.37-natIdt);
        out << "native sequence identity energy: " << natSeqWt<< endl;
        dt->updateNativeSequenceEnergy(natSeqWt);
    }



    vector<string> preSeqs;
    string preSeq = "";
    for(int i=1;i<=n;i++)
    {
        char x[10];
        char s[1000];
        sprintf(x,"000%d",i);
        string xx(x);
        xx = xx.substr(xx.length()-3,3);
        string fileName = output;
        if(n > 1)
            fileName = fileName.substr(0, output.length()-4) + "-" + xx + ".pdb";
        double score;
        if(preSeqs.size() >0 && cmd.specifiedOption("-nat")) {
            double preIdt = sequenceIdentity(natSeq, preSeq);
            if(abs(preIdt-natIdt) > 0.02) {
                natSeqWt = 5.0*(preIdt - natIdt);
                dt->updateNativeSequenceEnergy(natSeqWt);
            }
        }

        if(cmd.specifiedOption("-div")) {
            if(cmd.specifiedOption("-aaprop"))
                score = mc.mcRunWithIdtAACountRestrain(&unit,preSeqs,divCutoff,sigmaDiv,aacmc);
            else
                score = mc.mcRunWithIdtRestrain(&unit,preSeqs,divCutoff,sigmaDiv);
        } else {
            if(cmd.specifiedOption("-aaprop"))
                score = mc.mcRunWithAACountRestrain(&unit,aacmc);
            else
                score = mc.mcRun(&unit);
        }


        //RotSequence unit1(dt->rotGroups, aacr);
        //unit.copyValue(&unit1);

        double Es1, Erot, Eref;
        unit.calcEnergyComponents(dt, &Es1, &Erot, &Eref);

        string seq = unit.toAASequence();
        preSeqs.push_back(seq);
        preSeq = seq;
        sprintf(s,"%s %10.3f %10.3f",seq.c_str(), score, score-Eref);
        string result(s);
        out << result << endl;
        mc.printPDB(&unit, fileName);
    }

    out.close();
    delete dt;

    clock_t end = clock();
    cout << "time: " << (float)(end-start)/CLOCKS_PER_SEC << endl;

}


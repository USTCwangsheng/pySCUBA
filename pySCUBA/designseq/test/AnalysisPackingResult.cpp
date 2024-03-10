/*
 * AnalysisPackingResult.cpp
 *
 *  Created on: 2018Äê4ÔÂ24ÈÕ
 *      Author: notxp
 */


#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <time.h>
#include "designseq/ProteinRep.h"
#include "designseq/StructureInfo.h"
#include "geometry/xyz.h"
#include "designseq/AtomLib.h"
#include "designseq/StructureInfo.h"
#include "designseq/PackingAnalysis.h"
#include "designseq/AAProbabilityArray.h"
#include "designseq/CmdArgs.h"


using namespace std;
using namespace NSPdesignseq;
using namespace NSPgeometry;




int newAtom = 0;
int newRes = 0;
int newChain = 0;
int delAtom = 0;
int delRes = 0;
int delChain = 0;

clock_t START = clock();






void statPDB(PDB* pdb, AtomLib* atLib)
{
    vector<Residue*> resList = pdb->getResList();
    vector<PolarAtom*> polarAtomList;
    vector<int> resIDList;

    int resNum = resList.size();
    int i,j;
    Residue *res, *resPre;

    if(resList.size() > 0)
    {
        res = resList.at(0);
        if(res->hasAtom("O") && res->hasAtom("C"))
        {
            polarAtomList.push_back(new PolarAtom(res, "O"));
            resIDList.push_back(0);
        }
    }

    for(i=1;i<resNum;i++)
    {
        res = resList.at(i);
        resPre = resList.at(i-1);

        if(res->hasAtom("O") && res->hasAtom("C"))
        {
            polarAtomList.push_back(new PolarAtom(res, "O"));
            resIDList.push_back(i);
        }

        if(res->triName == "PRO")
            continue;

        if(res->hasAtom("N") && res->hasAtom("CA") && resPre->hasAtom("C"))
        {
            if(res->getAtom("N")->distance(*(resPre->getAtom("C"))) > 2)
                continue;
            polarAtomList.push_back(new PolarAtom(res,"N",resPre->getAtom("C")));
            resIDList.push_back(i);
        }
    }

    string resType;
    for(i=0;i<resNum;i++)
    {
        res = resList.at(i);
        resType = res->triName;
        if(resType == "ALA" || resType == "CYS" || resType == "PHE" ||
           resType == "GLY" || resType == "ILE" || resType == "LEU" ||
           resType == "MET" || resType == "PRO" || resType == "VAL")
            continue;
        if(!res->sidechainComplete(atLib))
            continue;
        if(resType == "ASP")
        {
            polarAtomList.push_back(new PolarAtom(res,"OD1"));
            polarAtomList.push_back(new PolarAtom(res,"OD2"));
            resIDList.push_back(i);
            resIDList.push_back(i);
        }
        else if(resType == "GLU")
        {
            polarAtomList.push_back(new PolarAtom(res,"OE1"));
            polarAtomList.push_back(new PolarAtom(res,"OE2"));
            resIDList.push_back(i);
            resIDList.push_back(i);
        }
        else if(resType == "HIS")
        {
            polarAtomList.push_back(new PolarAtom(res,"ND1"));
            polarAtomList.push_back(new PolarAtom(res,"NE2"));
            resIDList.push_back(i);
            resIDList.push_back(i);
        }
        else if(resType == "LYS")
        {
            polarAtomList.push_back(new PolarAtom(res,"NZ"));
            resIDList.push_back(i);
        }
        else if(resType == "ASN")
        {
            polarAtomList.push_back(new PolarAtom(res,"OD1"));
            polarAtomList.push_back(new PolarAtom(res,"ND2"));
            resIDList.push_back(i);
            resIDList.push_back(i);
        }
        else if(resType == "GLN")
        {
            polarAtomList.push_back(new PolarAtom(res,"OE1"));
            polarAtomList.push_back(new PolarAtom(res,"NE2"));
            resIDList.push_back(i);
            resIDList.push_back(i);
        }
        else if(resType == "ARG")
        {
            polarAtomList.push_back(new PolarAtom(res,"NE"));
            polarAtomList.push_back(new PolarAtom(res,"NH1"));
            polarAtomList.push_back(new PolarAtom(res,"NH2"));
            resIDList.push_back(i);
            resIDList.push_back(i);
            resIDList.push_back(i);
        }
        else if(resType == "SER")
        {
            polarAtomList.push_back(new PolarAtom(res,"OG"));
            resIDList.push_back(i);
        }
        else if(resType == "THR")
        {
            polarAtomList.push_back(new PolarAtom(res,"OG1"));
            resIDList.push_back(i);
        }
        else if(resType == "TRP")
        {
            polarAtomList.push_back(new PolarAtom(res,"NE1"));
            resIDList.push_back(i);
        }
        else if(resType == "TYR")
        {
            polarAtomList.push_back(new PolarAtom(res,"OH"));
            resIDList.push_back(i);
        }
    }
}

void statPDB2(PDB* pdb, AtomLib* atLib)
{
    vector<Residue*> resList = pdb->getResList();

    vector<PolarAtom*> polarAtomList;

    int resNum = resList.size();
    int i,j;
    Residue *res, *resPre;

    res = resList.at(0);
    if(res->hasAtom("O") && res->hasAtom("C"))
        polarAtomList.push_back(new PolarAtom(res, "O"));
    for(i=1;i<resNum;i++)
    {
        res = resList.at(i);
        resPre = resList.at(i-1);

        if(res->hasAtom("O") && res->hasAtom("C"))
            polarAtomList.push_back(new PolarAtom(res, "O"));

        if(res->triName == "PRO")
            continue;

        if(res->hasAtom("N") && res->hasAtom("CA") && resPre->hasAtom("C"))
        {
            if(res->getAtom("N")->distance(*(resPre->getAtom("C"))) > 2)
                continue;
            polarAtomList.push_back(new PolarAtom(res,"N",resPre->getAtom("C")));
        }
    }

    string resType;
    for(i=0;i<resNum;i++)
    {
        res = resList.at(i);
        resType = res->triName;
        if(resType == "ALA" || resType == "CYS" || resType == "PHE" ||
           resType == "GLY" || resType == "ILE" || resType == "LEU" ||
           resType == "MET" || resType == "PRO" || resType == "VAL")
            continue;
        if(!res->sidechainComplete(atLib))
            continue;
        if(resType == "ASP")
        {
            polarAtomList.push_back(new PolarAtom(res,"OD1"));
            polarAtomList.push_back(new PolarAtom(res,"OD2"));
        }
        else if(resType == "GLU")
        {
            polarAtomList.push_back(new PolarAtom(res,"OE1"));
            polarAtomList.push_back(new PolarAtom(res,"OE2"));
        }
        else if(resType == "HIS")
        {
            polarAtomList.push_back(new PolarAtom(res,"ND1"));
            polarAtomList.push_back(new PolarAtom(res,"NE2"));
        }
        else if(resType == "LYS")
        {
            polarAtomList.push_back(new PolarAtom(res,"NZ"));
        }
        else if(resType == "ASN")
        {
            polarAtomList.push_back(new PolarAtom(res,"OD1"));
            polarAtomList.push_back(new PolarAtom(res,"ND2"));
        }
        else if(resType == "GLN")
        {
            polarAtomList.push_back(new PolarAtom(res,"OE1"));
            polarAtomList.push_back(new PolarAtom(res,"NE2"));
        }
        else if(resType == "ARG")
        {
            polarAtomList.push_back(new PolarAtom(res,"NE"));
            polarAtomList.push_back(new PolarAtom(res,"NH1"));
            polarAtomList.push_back(new PolarAtom(res,"NH2"));
        }
        else if(resType == "SER")
        {
            polarAtomList.push_back(new PolarAtom(res,"OG"));
        }
        else if(resType == "THR")
        {
            polarAtomList.push_back(new PolarAtom(res,"OG1"));
        }
        else if(resType == "TRP")
        {
            polarAtomList.push_back(new PolarAtom(res,"NE1"));
        }
        else if(resType == "TYR")
        {
            polarAtomList.push_back(new PolarAtom(res,"OH"));
        }
    }
}

void statAll()
{
    ifstream input;
    string pdbListFile = "/home/xiongpeng/lib/pdb5/pdbList";

    input.open(pdbListFile.c_str(),ios::in);

    if (! input.is_open())
    {
        cout << "fail to open file " << pdbListFile << endl;
        exit (1);
    }
    string s, packOutput;

    AtomLib* atLib = new AtomLib();

    while(getline(input,s))
    {
        string third = s.substr(2,1);
        PDB* pdb = new PDB("/home/xiongpeng/lib/pdb5/" + third + "/"+s+".pdb", s);
        statPDB(pdb,atLib);
        delete pdb;
    }


}



/*

void anaHBondList(string listTag)
{
    ResName rn;

    ifstream input;
    string pdbListFile = "/export/home/s2982206/pack/dataset/list1000-"+listTag;

    input.open(pdbListFile.c_str(),ios::in);

    if (! input.is_open())
    {
        cout << "fail to open file " << pdbListFile << endl;
        exit (1);
    }
    vector<string> pdbList;
    string s, packOutput;
    while(getline(input,s))
    {
        pdbList.push_back(s);
    }

    string output = "/export/home/s2982206/pack/result/list1000-"+listTag+".hbstat";
    FILE* pFile;
    pFile = fopen(output.c_str(), "w");



    AtomLib* atLib = new AtomLib();

    set<string> setN;
    set<string> setP;

    int countNat[54];
    int countPack[54];
    for(int i=0;i<54;i++)
    {
        countNat[i] = 0;
        countPack[i] = 0;
    }

    for(size_t i=0;i<pdbList.size();i++)
    {

        s = pdbList.at(i);
        cout << s << endl;
        PDB* pdbA = new PDB("/export/home/s2982206/pack/dataset/pdb/"+s+".pdb", s);
        packOutput = "/export/home/s2982206/pack/dataset/pack/"+s+"-x.pdb";


        PDB* pdbB = new PDB(packOutput, s);

        hydrogenBondCount(pdbA,rsp,atLib,countNat, setN);
        cout << "finish nat: " << endl;
        hydrogenBondCount(pdbB,rsp,atLib,countPack, setP);
        delete pdbA;
        delete pdbB;
    }

    set<string>::iterator it,it2;


    int inNat = setN.size();
    int inPack = setP.size();
    int inBoth = 0;

    for(it = setN.begin();it!=setN.end();it++)
    {
        string tag = *it;
        it2 = setP.find(tag);
        if(it2 != setP.end())
            inBoth++;
    }

    string doners[] = {"BB-N", "HIS-NDE", "ARG-NH", "ARG-NE", "LYS-NZ", "TRP-NE", "TS-OG", "TYR-OH", "NQ-NDE"};
    string acceptors[] = {"BB-O", "DE-O", "NQ-O", "TS-OG", "TYR-OH", "HIS-NDE" };

    for(int i=0;i<54;i++)
    {
        int idAcc = i%6;
        int idDoner = i/6;
        int pairType = idDoner*6+idAcc;
        if(pairType != i)
        {
            cout << "wrong answer" << endl;
        }
        fprintf(pFile,"%-7s %-7s countNat: %7d countPack: %7d\n", doners[idDoner].c_str(), acceptors[idAcc].c_str(), countNat[pairType], countPack[pairType] );
    }

    fprintf(pFile, "inNat: %6d inPack: %6d inBoth: %6d\n",inNat, inPack, inBoth);
    fclose(pFile);
}

void anaHBond(string paraTag)
{
    ResName rn;

    ifstream input;
    string pdbListFile = "/export/home/s2982206/pack/dataset/list1000";

    input.open(pdbListFile.c_str(),ios::in);

    if (! input.is_open())
    {
        cout << "fail to open file " << pdbListFile << endl;
        exit (1);
    }
    vector<string> pdbList;
    string s, packOutput;
    while(getline(input,s))
    {
        pdbList.push_back(s);
    }

    string output = "/export/home/s2982206/pack/result/"+paraTag+".hbstat";
    FILE* pFile;
    pFile = fopen(output.c_str(), "w");

    vector<float> distanceList;
    vector<float> saiList;
    vector<string> distTypeList;

    vector<float> chi1DiffList;
    vector<float> chi1DiffSaiList;
    vector<string> chi1DiffTypeList;

    vector<string> inBoth;
    vector<string> onlyA;
    vector<string> onlyB;

    AtomLib* atLib = new AtomLib();
    ResSasaPoints* rsp = new ResSasaPoints();

    int countNat[54];
    int countPack[54];
    for(int i=0;i<54;i++)
    {
        countNat[i] = 0;
        countPack[i] = 0;
    }

    set<string> setN;
    set<string> setP;

    for(size_t i=0;i<pdbList.size();i++)
    {

        s = pdbList.at(i);
        cout << s << endl;
        PDB* pdbA = new PDB("/export/home/s2982206/pack/dataset/pdb/"+s+".pdb", s);
        packOutput = "/export/home/s2982206/pack/dataset/pack/"+s+"-" + paraTag + ".pdb";


        PDB* pdbB = new PDB(packOutput, s);

        hydrogenBondCount(pdbA,rsp,atLib,countNat, setN);
        hydrogenBondCount(pdbB,rsp,atLib,countPack, setP);
        delete pdbA;
        delete pdbB;
    }

    string doners[] = {"BB-N", "HIS-NDE", "ARG-NH", "ARG-NE", "LYS-NZ", "TRP-NE", "TS-OG", "TYR-OH", "NQ-NDE"};
    string acceptors[] = {"BB-O", "DE-O", "NQ-O", "TS-OG", "TYR-OH", "HIS-NDE" };

    for(int i=0;i<54;i++)
    {
        int idDoner = i%9;
        int idAcc = i/9;
        int pairType = idDoner*6+idAcc;
        printf("%-6s %-6s countNat: %7d countPack: %7d\n", doners[idDoner].c_str(), acceptors[idAcc].c_str(), countNat[pairType], countPack[pairType] );
    }
}
*/
/*
void anaPackingResult(string paraTag, string triName)
{
    ResName rn;

    ifstream input;
    string pdbListFile = "/export/home/s2982206/pack/dataset/list100";

    input.open(pdbListFile.c_str(),ios::in);

    if (! input.is_open())
    {
        cout << "fail to open file " << pdbListFile << endl;
        exit (1);
    }
    vector<string> pdbList;
    string s, packOutput;
    while(getline(input,s))
    {
        pdbList.push_back(s);
    }

    vector<float> distanceList;

    vector<float> chi1DiffList;

    AtomLib* atLib = new AtomLib();

    for(size_t i=0;i<pdbList.size();i++)
    {

        s = pdbList.at(i);
        //cout << s << endl;
        //PDB* pdbA = new PDB("/export/home/s2982206/pack/dataset/pdb/"+s+".pdb", s);

        PDB* pdbA = new PDB("/export/home/s2982206/pack/dataset/pack/"+s+"-x.pdb", s);
        packOutput = "/export/home/s2982206/pack/dataset/pack/"+s+"-" + paraTag + ".pdb";


        PDB* pdbB = new PDB(packOutput, s);

        atomicRMS(pdbA,pdbB,atLib,distanceList,triName);

        chi1Diff(pdbA,pdbB,atLib,chi1DiffList, triName);

        delete pdbA;
        delete pdbB;
    }


    vector<float> aaDistList[20];


    float averageRMS = squareAverage(distanceList);


    int chi1Succ = 0;
    int chi1Tot = 0;


    for(size_t i=0;i<chi1DiffList.size();i++)
    {
        chi1Tot++;
        if(chi1DiffList.at(i) < 40)
        {
            chi1Succ ++;
        }
    }

    printf("rms: %5.3f chi1: %6.4f\n",averageRMS, chi1Succ*1.0/chi1Tot);

}
*/

void anaPDB(string natPDBFile, string packPDBFile, string output){
    vector<float> distanceList;
    vector<float> saiList;
    vector<string> distTypeList;

    vector<float> chi1DiffList;
    vector<float> chi1DiffSaiList;
    vector<string> chi1DiffTypeList;

    ofstream out(output, ios::out);

    PDB* pdbA = new PDB(natPDBFile, "unk");
    PDB* pdbB = new PDB(packPDBFile, "unk");

    AtomLib* atLib = new AtomLib();
    atomicRMS(pdbA,pdbB,atLib,distanceList,saiList, distTypeList);
    chi1Diff(pdbA,pdbB,atLib, chi1DiffList, chi1DiffSaiList, chi1DiffTypeList);




    vector<string> inBoth;
    vector<string> onlyA;
    vector<string> onlyB;
    disulfideBondDiff(pdbA,pdbB,inBoth,onlyA,onlyB);
    delete pdbA;
    delete pdbB;


    char s[100];
    for(int i=0;i<distanceList.size();i++){
        sprintf(s,"dist: %s %5.3f %6.4f\n",distTypeList[i].c_str(), saiList[i], distanceList[i]);
        string x(s);
        out << x ;
    }

    for(int i=0;i<chi1DiffList.size();i++){
        sprintf(s,"chi1: %s %5.3f %8.3f\n",chi1DiffTypeList[i].c_str(), chi1DiffSaiList[i], chi1DiffList[i]);
        string x(s);
        out << x;
    }

    out << "DS " << inBoth.size() << " " << onlyA.size() << " " << onlyB.size() << "\n";

    out.close();

}

void anaPackingResult(string paraTag, string output)
{
    ResName rn;

    ifstream input;
    string pdbListFile = "/home/xiongpeng/packing/natPDB/list/training";

    input.open(pdbListFile.c_str(),ios::in);

    if (! input.is_open())
    {
        cout << "fail to open file " << pdbListFile << endl;
        exit (1);
    }
    vector<string> pdbList;
    string s, packOutput;
    while(getline(input,s))
    {
        pdbList.push_back(s);
    }


    FILE* pFile;
    pFile = fopen(output.c_str(), "w");

    vector<float> distanceList;
    vector<float> saiList;
    vector<string> distTypeList;

    vector<float> chi1DiffList;
    vector<float> chi1DiffSaiList;
    vector<string> chi1DiffTypeList;

    vector<string> inBoth;
    vector<string> onlyA;
    vector<string> onlyB;

    AtomLib* atLib = new AtomLib();

    for(size_t i=0;i<pdbList.size();i++)
    {

        s = pdbList.at(i);

        string third = s.substr(2,1);

        PDB* pdbA = new PDB("/home/xiongpeng/packing/natPDB/flip/"+s+".pdb", s);
        packOutput = "/home/xiongpeng/packing/packing/"+s+"-" + paraTag + ".pdb";


        PDB* pdbB = new PDB(packOutput, s);


        atomicRMS(pdbA,pdbB,atLib,distanceList,saiList, distTypeList);

        if(distTypeList.size() != distTypeList.size())
        {
            cerr << s << " size error" << endl;
            exit(0);
        }


        chi1Diff(pdbA,pdbB,atLib, chi1DiffList, chi1DiffSaiList, chi1DiffTypeList);
        if(chi1DiffList.size() != chi1DiffTypeList.size())
        {
            cerr << s << " chi size error " << endl;
            exit(0);
        }



        disulfideBondDiff(pdbA,pdbB,inBoth,onlyA,onlyB);
        delete pdbA;
        delete pdbB;
    }


    vector<float> aaDistList[20];

    vector<float> coreDistList;
    for(size_t i=0;i<distanceList.size();i++)
    {
        int type = rn.triToInt(distTypeList.at(i));
        aaDistList[type].push_back(distanceList.at(i));
        if(saiList.at(i) < 0.33)
            coreDistList.push_back(distanceList.at(i));
    }
    float averageRMS = squareAverage(distanceList);
    float averageCoreRMS = squareAverage(coreDistList);


    int chi1Succ = 0;
    int chi1Tot = 0;
    int chi1SuccAA[20] = {0};
    int chi1AATot[20] = {0};

    int chi1CoreSucc = 0;
    int chi1CoreTot = 0;
    int chi1CoreSuccAA[20] = {0};
    int chi1CoreAATot[20] = {0};

    for(size_t i=0;i<chi1DiffList.size();i++)
    {
        int type = rn.triToInt(chi1DiffTypeList.at(i));
        chi1Tot++;
        chi1AATot[type] ++;
        if(chi1DiffList.at(i) < 40)
        {
            chi1Succ ++;
            chi1SuccAA[type] ++;
        }

        if(chi1DiffSaiList.at(i) < 0.33)
        {
            chi1CoreTot++;
            chi1CoreAATot[type] ++;
            if(chi1DiffList.at(i) < 40)
            {
                chi1CoreSucc ++;
                chi1CoreSuccAA[type] ++;
            }
        }
    }

    cout << paraTag << endl;
    fprintf(pFile,"TOT RMS: %5.3f  CoreRMS: %5.3f dsTot: %3d dsNat: %3d dsPack: %3d\n",averageRMS, averageCoreRMS, inBoth.size(), onlyA.size(), onlyB.size());

    fprintf(pFile,"TOT chi1 succecc rate: %6.4f Core success rate: %6.4f\n", chi1Succ*1.0/chi1Tot , chi1CoreSucc*1.0/chi1CoreTot);
    for(int i=0;i<20;i++)
    {
        if(chi1AATot[i] > 0)
        {
            string type = rn.intToTri(i);
        //    cout << aaDistList[i].size() << endl;
        //    cout << aaDistList[i].at(0) << " " << aaDistList[i].at(1) << endl;
            float rms = squareAverage(aaDistList[i]);
            fprintf(pFile,"%s %6.4f %6.4f\n",type.c_str(),chi1SuccAA[i]*1.0/chi1AATot[i], rms);
        }
    }

    fclose(pFile);
}


int main(int argc, char** argv)
{

    anaPackingResult(argv[1], argv[2]);


}


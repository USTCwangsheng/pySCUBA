/*
 * PackingAnalysis.cpp
 *
 *  Created on: 2018Äê4ÔÂ24ÈÕ
 *      Author: notxp
 */

#include "designseq/PackingAnalysis.h"

namespace NSPdesignseq {

float squareAverage(vector<float>& list){
    if(list.size()==0)
        return 0;
    double tot = 0;
    for(int i=0;i<list.size();i++){
        tot += list[i]*list[i];

    }
    return sqrt(tot/list.size());
}

vector<float> distListLeu(Residue* resA, Residue* resB)
{
    if(resA->triName != "LEU")
    {
        cerr << " only for LEU: " << resA->triName << endl;
    }
    double e1 = 0;
    double e2 = 0;
    vector<float> list1;
    vector<float> list2;

    list1.push_back(resA->getAtom("CB")->distance(* (resB->getAtom("CB"))));
    list1.push_back(resA->getAtom("CG")->distance(* (resB->getAtom("CG"))));
    list1.push_back(resA->getAtom("CD1")->distance(* (resB->getAtom("CD1"))));
    list1.push_back(resA->getAtom("CD2")->distance(* (resB->getAtom("CD2"))));

    list2.push_back(resA->getAtom("CB")->distance(* (resB->getAtom("CB"))));
    list2.push_back(resA->getAtom("CG")->distance(* (resB->getAtom("CG"))));
    list2.push_back(resA->getAtom("CD1")->distance(* (resB->getAtom("CD2"))));
    list2.push_back(resA->getAtom("CD2")->distance(* (resB->getAtom("CD1"))));

    e1 = squareAverage(list1);
    e2 = squareAverage(list2);


    if(e1 < e2)
        return list1;
    else
        return list2;

}

vector<float> distListVal(Residue* resA, Residue* resB)
{
    if(resA->triName != "VAL")
    {
        cerr << " only for VAL: " << resA->triName << endl;
    }
    double e1 = 0;
    double e2 = 0;
    vector<float> list1;
    vector<float> list2;

    list1.push_back(resA->getAtom("CB")->distance(* (resB->getAtom("CB"))));
    list1.push_back(resA->getAtom("CG1")->distance(* (resB->getAtom("CG1"))));
    list1.push_back(resA->getAtom("CG2")->distance(* (resB->getAtom("CG2"))));

    list2.push_back(resA->getAtom("CB")->distance(* (resB->getAtom("CB"))));
    list2.push_back(resA->getAtom("CG1")->distance(* (resB->getAtom("CG2"))));
    list2.push_back(resA->getAtom("CG2")->distance(* (resB->getAtom("CG1"))));

    e1 = squareAverage(list1);
    e2 = squareAverage(list2);


    if(e1 < e2)
        return list1;
    else
        return list2;
}

vector<float> distListPhe(Residue* resA, Residue* resB)
{
    if(resA->triName != "PHE")
    {
        cerr << " only for PHE: " << resA->triName << endl;
    }
    double e1 = 0;
    double e2 = 0;
    vector<float> list1;
    vector<float> list2;

    list1.push_back(resA->getAtom("CB")->distance(* (resB->getAtom("CB"))));
    list1.push_back(resA->getAtom("CG")->distance(* (resB->getAtom("CG"))));
    list1.push_back(resA->getAtom("CD1")->distance(* (resB->getAtom("CD1"))));
    list1.push_back(resA->getAtom("CD2")->distance(* (resB->getAtom("CD2"))));
    list1.push_back(resA->getAtom("CE1")->distance(* (resB->getAtom("CE1"))));
    list1.push_back(resA->getAtom("CE2")->distance(* (resB->getAtom("CE2"))));
    list1.push_back(resA->getAtom("CZ")->distance(* (resB->getAtom("CZ"))));

    list2.push_back(resA->getAtom("CB")->distance(* (resB->getAtom("CB"))));
    list2.push_back(resA->getAtom("CG")->distance(* (resB->getAtom("CG"))));
    list2.push_back(resA->getAtom("CD1")->distance(* (resB->getAtom("CD2"))));
    list2.push_back(resA->getAtom("CD2")->distance(* (resB->getAtom("CD1"))));
    list2.push_back(resA->getAtom("CE1")->distance(* (resB->getAtom("CE2"))));
    list2.push_back(resA->getAtom("CE2")->distance(* (resB->getAtom("CE1"))));
    list2.push_back(resA->getAtom("CZ")->distance(* (resB->getAtom("CZ"))));

    e1 = squareAverage(list1);
    e2 = squareAverage(list2);


    if(e1 < e2)
        return list1;
    else
        return list2;
}

vector<float> distListTyr(Residue* resA, Residue* resB)
{
    if(resA->triName != "TYR")
    {
        cerr << " only for TYR: " << resA->triName << endl;
    }
    double e1 = 0;
    double e2 = 0;
    vector<float> list1;
    vector<float> list2;

    list1.push_back(resA->getAtom("CB")->distance(* (resB->getAtom("CB"))));
    list1.push_back(resA->getAtom("CG")->distance(* (resB->getAtom("CG"))));
    list1.push_back(resA->getAtom("CD1")->distance(* (resB->getAtom("CD1"))));
    list1.push_back(resA->getAtom("CD2")->distance(* (resB->getAtom("CD2"))));
    list1.push_back(resA->getAtom("CE1")->distance(* (resB->getAtom("CE1"))));
    list1.push_back(resA->getAtom("CE2")->distance(* (resB->getAtom("CE2"))));
    list1.push_back(resA->getAtom("CZ")->distance(* (resB->getAtom("CZ"))));
    list1.push_back(resA->getAtom("OH")->distance(* (resB->getAtom("OH"))));

    list2.push_back(resA->getAtom("CB")->distance(* (resB->getAtom("CB"))));
    list2.push_back(resA->getAtom("CG")->distance(* (resB->getAtom("CG"))));
    list2.push_back(resA->getAtom("CD1")->distance(* (resB->getAtom("CD2"))));
    list2.push_back(resA->getAtom("CD2")->distance(* (resB->getAtom("CD1"))));
    list2.push_back(resA->getAtom("CE1")->distance(* (resB->getAtom("CE2"))));
    list2.push_back(resA->getAtom("CE2")->distance(* (resB->getAtom("CE1"))));
    list2.push_back(resA->getAtom("CZ")->distance(* (resB->getAtom("CZ"))));
    list2.push_back(resA->getAtom("OH")->distance(* (resB->getAtom("OH"))));

    e1 = squareAverage(list1);
    e2 = squareAverage(list2);


    if(e1 < e2)
        return list1;
    else
        return list2;
}

vector<float> distListAsp(Residue* resA, Residue* resB)
{
    if(resA->triName != "ASP")
    {
        cerr << " only for ASP: " << resA->triName << endl;
    }
    double e1 = 0;
    double e2 = 0;
    vector<float> list1;
    vector<float> list2;

    list1.push_back(resA->getAtom("CB")->distance(* (resB->getAtom("CB"))));
    list1.push_back(resA->getAtom("CG")->distance(* (resB->getAtom("CG"))));
    list1.push_back(resA->getAtom("OD1")->distance(* (resB->getAtom("OD1"))));
    list1.push_back(resA->getAtom("OD2")->distance(* (resB->getAtom("OD2"))));

    list2.push_back(resA->getAtom("CB")->distance(* (resB->getAtom("CB"))));
    list2.push_back(resA->getAtom("CG")->distance(* (resB->getAtom("CG"))));
    list2.push_back(resA->getAtom("OD1")->distance(* (resB->getAtom("OD2"))));
    list2.push_back(resA->getAtom("OD2")->distance(* (resB->getAtom("OD1"))));

    e1 = squareAverage(list1);
    e2 = squareAverage(list2);


    if(e1 < e2)
        return list1;
    else
        return list2;
}

vector<float> distListGlu(Residue* resA, Residue* resB)
{
    if(resA->triName != "GLU")
    {
        cerr << " only for GLU: " << resA->triName << endl;
    }
    double e1 = 0;
    double e2 = 0;
    vector<float> list1;
    vector<float> list2;



    list1.push_back(resA->getAtom("CB")->distance(* (resB->getAtom("CB"))));
    list1.push_back(resA->getAtom("CG")->distance(* (resB->getAtom("CG"))));
    list1.push_back(resA->getAtom("CD")->distance(* (resB->getAtom("CD"))));
    list1.push_back(resA->getAtom("OE1")->distance(* (resB->getAtom("OE1"))));
    list1.push_back(resA->getAtom("OE2")->distance(* (resB->getAtom("OE2"))));

    list2.push_back(resA->getAtom("CB")->distance(* (resB->getAtom("CB"))));
    list2.push_back(resA->getAtom("CG")->distance(* (resB->getAtom("CG"))));
    list2.push_back(resA->getAtom("CD")->distance(* (resB->getAtom("CD"))));
    list2.push_back(resA->getAtom("OE1")->distance(* (resB->getAtom("OE2"))));
    list2.push_back(resA->getAtom("OE2")->distance(* (resB->getAtom("OE1"))));

    e1 = squareAverage(list1);
    e2 = squareAverage(list2);


    if(e1 < e2)
        return list1;
    else
        return list2;
}


int atomicRMS(PDB* pdbA, PDB* pdbB, AtomLib* atLib, vector<float>& distanceList, string triName)
{
    /*
     * return total size chain atom number
     */

    int n = 0;
    vector<Residue*> tmpAList = pdbA->getValidResList();
    vector<Residue*> tmpBList = pdbB->getValidResList();

    if(tmpAList.size() != tmpBList.size())
    {
        cout << pdbA->getPDBID() << " resList not equal" << endl;
    }

    vector<Residue*> resAList;
    vector<Residue*> resBList;

    for(size_t i=0;i<tmpAList.size();i++)
    {
        if(tmpAList.at(i)->resID != tmpBList.at(i)->resID)
        {
            cout << pdbA->getPDBID() << endl;
            cout << tmpAList.at(i)->resID << " " << tmpBList.at(i)->resID << endl;
        }
        Residue* resA = tmpAList.at(i);
        Residue* resB = tmpBList.at(i);
        if(!resA->hasThreeCoreAtoms()) continue;
        if(!resB->hasThreeCoreAtoms()) continue;
        if(!resA->hasAtom("O")) continue;
        if(!resB->hasAtom("O")) continue;

        if(resA->altLoc != ' ') continue;
        if(resB->altLoc != ' ') continue;
        resAList.push_back(resA);
        resBList.push_back(resB);
    }

    if(resAList.size() != resBList.size())
    {
        cerr << "residue list size not equal: " << resAList.size() << " " << resBList.size() << endl;
        exit(1);
    }

    for(size_t i=0;i<resAList.size();i++)
    {
        Residue* resA = resAList.at(i);
        Residue* resB = resBList.at(i);
        if(resA->triName != resB->triName)
        {
            cerr << "PDBID: " << pdbA->getPDBID() << " Residue type not consistent: " << resA->resID << " " << resA->triName << " " << resB->triName << endl;
            exit(1);
        }


        if(resA->triName == triName)
        {
            vector<string>* scAtomNames = atLib->getAminoAcidSidechainAtomNames(resA->triName);
            string atomName;

            for(size_t j=0;j<scAtomNames->size();j++)
            {
                atomName = scAtomNames->at(j);
                if(!resA->hasAtom(atomName) || !resB->hasAtom(atomName))
                {
                    cerr << "sidechain not complete " << pdbA->getPDBID() << "  res: " << resA->resID << " type: " << resA->triName << " " << atomName << endl;
                    continue;
                }
                else
                {
                    n++;
                    distanceList.push_back(resA->getAtom(atomName)->distance(*resB->getAtom(atomName)));
                }
            }
        }
    }
    return n;
}

int atomicRMS(PDB* pdbA, PDB* pdbB, AtomLib* atLib, vector<float>& distanceList, vector<float>& saiList, vector<string>& typeList)
{
    /*
     * return total size chain atom number
     */


    vector<Residue*> tmpAList = pdbA->getValidResList();
    vector<Residue*> tmpBList = pdbB->getValidResList();

    if(tmpAList.size() != tmpBList.size())
    {
        cout << "list size not equal: " << pdbA->getPDBID() << " " << tmpAList.size() << " " << tmpBList.size() << endl;
        return 0;
    }

    int num = 0;
    vector<Residue*> resAList;
    vector<Residue*> resBList;

    for(size_t i=0;i<tmpAList.size();i++)
    {
        Residue* resA = tmpAList.at(i);
        Residue* resB = tmpBList.at(i);
        if(!resA->hasThreeCoreAtoms()) continue;
        if(!resB->hasThreeCoreAtoms()) continue;
        if(!resA->hasAtom("O")) continue;
        if(!resB->hasAtom("O")) continue;
        if(resA->altLoc != ' ') continue;
        if(resB->altLoc != ' ') continue;

        resAList.push_back(resA);
        resBList.push_back(resB);
    }


    StructureInfo si(pdbB);
    SasaPSD rsp;
    si.updateSAI(&rsp);
    vector<float> pdbSaiList = si.getSaiList();
    if(resAList.size() != resBList.size())
    {
        cerr << "residue list size not equal: " << resAList.size() << " " << resBList.size() << endl;
        exit(1);
    }

    if(resAList.size()!= pdbSaiList.size())
    {
        cout << pdbA->getPDBID() << " sai size error :" << resAList.size() << " " << pdbSaiList.size() << endl;
    }

    for(size_t i=0;i<resAList.size();i++)
    {

        Residue* resA = resAList.at(i);
        Residue* resB = resBList.at(i);
        if(resA->altLoc != ' ') continue;
        if(resB->altLoc != ' ') continue;
        if(!resA->sidechainComplete(atLib) || !resB->sidechainComplete(atLib))
        {
            cerr << "sidechain not complete :" << pdbA->getPDBID() << " " << resA->getResID() << " " << resA->triName << endl;
            continue;
        }

        if(resA->triName != resB->triName)
        {
            cerr << "PDBID: " << pdbA->getPDBID() << " Residue type not consistent: " << resA->resID << " " << resA->triName << " " << resB->triName << endl;
            continue;
        }


        if(resA->triName == "LEU")
        {
            vector<float> list = distListLeu(resA,resB);
            for(size_t k=0;k<list.size();k++)
            {
                num++;
                distanceList.push_back(list.at(k));
                typeList.push_back("LEU");
                saiList.push_back(pdbSaiList.at(i));
            }
            continue;
        }
        else if(resA->triName == "VAL")
        {
            vector<float> list = distListVal(resA,resB);
            for(size_t k=0;k<list.size();k++)
            {
                num++;
                distanceList.push_back(list.at(k));
                typeList.push_back("VAL");
                saiList.push_back(pdbSaiList.at(i));
            }
            continue;
        }
        else if(resA->triName == "PHE")
        {
            vector<float> list = distListPhe(resA,resB);
            for(size_t k=0;k<list.size();k++)
            {
                num++;
                distanceList.push_back(list.at(k));
                typeList.push_back("PHE");
                saiList.push_back(pdbSaiList.at(i));
            }
            continue;
        }
        else if(resA->triName == "TYR")
        {
            vector<float> list = distListTyr(resA,resB);
            for(size_t k=0;k<list.size();k++)
            {
                num++;
                distanceList.push_back(list.at(k));
                typeList.push_back("TYR");
                saiList.push_back(pdbSaiList.at(i));
            }
            continue;
        }
        else if(resA->triName == "ASP")
        {
            vector<float> list = distListAsp(resA,resB);
            for(size_t k=0;k<list.size();k++)
            {
                num++;
                distanceList.push_back(list.at(k));
                typeList.push_back("ASP");
                saiList.push_back(pdbSaiList.at(i));
            }
            continue;
        }
        else if(resA->triName == "GLU")
        {

            vector<float> list = distListGlu(resA,resB);

            for(size_t k=0;k<list.size();k++)
            {
                num++;

                distanceList.push_back(list.at(k));

                typeList.push_back("GLU");

                saiList.push_back(pdbSaiList.at(i));

            }
            continue;
        }



        vector<string>* scAtomNames = atLib->getAminoAcidSidechainAtomNames(resA->triName);
        string atomName;

        for(size_t j=0;j<scAtomNames->size();j++)
        {
            atomName = scAtomNames->at(j);
            if(!resA->hasAtom(atomName) || !resB->hasAtom(atomName))
            {
                cerr << "sidechain not complete: " <<pdbA->getPDBID() << " " <<  resA->resID << " type: " << resA->triName << " " << atomName << endl;
                continue;
            }
            else
            {
                num++;
                distanceList.push_back(resA->getAtom(atomName)->distance(*resB->getAtom(atomName)));
                saiList.push_back(pdbSaiList.at(i));
                typeList.push_back(resA->triName);
            }
        }
    }
    return num;
}

void chi1Diff(PDB* pdbA, PDB* pdbB, AtomLib* atLib, vector<float>& chi1DiffList, string triName)
{
    ResName rn;
    vector<Residue*> tmpAList = pdbA->getValidResList();
    vector<Residue*> tmpBList = pdbB->getValidResList();
    if(tmpAList.size() != tmpBList.size())
    {
        cout << "list size not equal: " << pdbA->getPDBID() << endl;
        return;
    }

    int num = 0;
    vector<Residue*> resAList;
    vector<Residue*> resBList;

    for(size_t i=0;i<tmpAList.size();i++)
    {
        Residue* resA = tmpAList.at(i);
        Residue* resB = tmpBList.at(i);
        if(!resA->hasThreeCoreAtoms()) continue;
        if(!resB->hasThreeCoreAtoms()) continue;
        if(resA->altLoc != ' ') continue;
        if(resB->altLoc != ' ') continue;
        resAList.push_back(resA);
        resBList.push_back(resB);
    }


    if(resAList.size() != resBList.size())
    {
        cerr << "residue list size not equal: " << resAList.size() << " " << resBList.size() << endl;
        exit(1);
    }

    float chi1A, chi1B, chi1Diff;
    for(size_t i=0;i<resAList.size();i++)
    {
        Residue* resA = resAList.at(i);
        Residue* resB = resBList.at(i);
        if(resA->triName != resB->triName)
        {
            cerr << "Residue type not consistent: " << resA->resID << " " << resA->triName << " " << resB->triName << endl;
            exit(1);
        }
        if(resA->triName != triName)
            continue;

        if(!resA->sidechainComplete(atLib) || !resB->sidechainComplete(atLib))
            continue;
        if(!resA->hasThreeCoreAtoms() || !resB->hasThreeCoreAtoms())
            continue;
        if(!resA->hasAtom("O")) continue;
        if(!resB->hasAtom("O")) continue;
        if(resA->altLoc != ' ') continue;
        if(resB->altLoc != ' ') continue;

        int chiNum = rn.chiNum(resA->triName);
        if(chiNum == 0)
            continue;
        vector<string>* scAtomNames = atLib->getAminoAcidSidechainAtomNames(resA->triName);
        chi1A = dihedral(resA->getAtom("N")->getCoord(), resA->getAtom("CA")->getCoord(), resA->getAtom("CB")->getCoord(), resA->getAtom(scAtomNames->at(1))->getCoord());
        chi1B = dihedral(resB->getAtom("N")->getCoord(), resB->getAtom("CA")->getCoord(), resB->getAtom("CB")->getCoord(), resB->getAtom(scAtomNames->at(1))->getCoord());
        chi1Diff = chi1A - chi1B;
        if(chi1Diff > 180)
            chi1Diff -= 360;
        if(chi1Diff < -180)
            chi1Diff += 360;
        chi1DiffList.push_back(abs(chi1Diff));
    }

}


void chi1Diff(PDB* pdbA, PDB* pdbB, AtomLib* atLib,vector<float>& chi1DiffList, vector<float>& saiList, vector<string>& resTypeList)
{
    ResName rn;
    vector<Residue*> tmpAList = pdbA->getValidResList();
    vector<Residue*> tmpBList = pdbB->getValidResList();
    if(tmpAList.size() != tmpBList.size())
    {
        cout << "list size not equal: " << pdbA->getPDBID() << endl;
        return ;
    }

    int num = 0;
    vector<Residue*> resAList;
    vector<Residue*> resBList;

    for(size_t i=0;i<tmpAList.size();i++)
    {
        Residue* resA = tmpAList.at(i);
        Residue* resB = tmpBList.at(i);
        if(!resA->hasThreeCoreAtoms()) continue;
        if(!resB->hasThreeCoreAtoms()) continue;
        if(!resA->hasAtom("O")) continue;
        if(!resB->hasAtom("O")) continue;
        if(resA->altLoc != ' ') continue;
        if(resB->altLoc != ' ') continue;

        resAList.push_back(resA);
        resBList.push_back(resB);
    }

    StructureInfo si(pdbB);
    SasaPSD rsp;
    si.updateSAI(&rsp);
    vector<float> pdbSaiList = si.getSaiList();
    if(resAList.size() != resBList.size())
    {
        cerr << "residue list size not equal: " << resAList.size() << " " << resBList.size() << endl;
        exit(1);
    }

    float chi1A, chi1B, chi1Diff;

    for(size_t i=0;i<resAList.size();i++)
    {
        Residue* resA = resAList.at(i);
        Residue* resB = resBList.at(i);
        if(resA->altLoc != ' ') continue;
        if(resB->altLoc != ' ') continue;
        if(resA->triName != resB->triName)
        {
            cerr << "Residue type not consistent: " << resA->resID << " " << resA->triName << " " << resB->triName << endl;
            continue;
        }

        if(!resA->sidechainComplete(atLib) || !resB->sidechainComplete(atLib))
            continue;
        if(!resA->hasThreeCoreAtoms() || !resB->hasThreeCoreAtoms())
            continue;
        if(!resA->hasAtom("O")) continue;
        if(!resB->hasAtom("O")) continue;
        if(resA->altLoc != ' ') continue;
        if(resB->altLoc != ' ') continue;

        int chiNum = rn.chiNum(resA->triName);
        if(chiNum == 0)
            continue;
        vector<string>* scAtomNames = atLib->getAminoAcidSidechainAtomNames(resA->triName);
        chi1A = dihedral(resA->getAtom("N")->getCoord(), resA->getAtom("CA")->getCoord(), resA->getAtom("CB")->getCoord(), resA->getAtom(scAtomNames->at(1))->getCoord());
        chi1B = dihedral(resB->getAtom("N")->getCoord(), resB->getAtom("CA")->getCoord(), resB->getAtom("CB")->getCoord(), resB->getAtom(scAtomNames->at(1))->getCoord());
        chi1Diff = chi1A - chi1B;
        if(chi1Diff > 180)
            chi1Diff -= 360;
        if(chi1Diff < -180)
            chi1Diff += 360;
        chi1DiffList.push_back(abs(chi1Diff));
        saiList.push_back(pdbSaiList.at(i));
        resTypeList.push_back(resA->triName);
    }
}

void chi2Diff(PDB* pdbA, PDB* pdbB,AtomLib* atLib, vector<float>& chi2DiffList, vector<float>& saiList, vector<string>& resTypeList)
{
    ResName rn;
    vector<Residue*> tmpAList = pdbA->getValidResList();
    vector<Residue*> tmpBList = pdbB->getValidResList();
    if(tmpAList.size() != tmpBList.size())
    {
        cout << "list size not equal: " << pdbA->getPDBID() << endl;
        return ;
    }

    int num = 0;
    vector<Residue*> resAList;
    vector<Residue*> resBList;

    for(size_t i=0;i<tmpAList.size();i++)
    {
        Residue* resA = tmpAList.at(i);
        Residue* resB = tmpBList.at(i);
        if(!resA->hasThreeCoreAtoms()) continue;
        if(!resB->hasThreeCoreAtoms()) continue;
        if(!resA->hasAtom("O")) continue;
        if(!resB->hasAtom("O")) continue;
        if(resA->altLoc != ' ') continue;
        if(resB->altLoc != ' ') continue;

        resAList.push_back(resA);
        resBList.push_back(resB);
    }
    StructureInfo si(pdbB);
    SasaPSD rsp;
    si.updateSAI(&rsp);
    vector<float> pdbSaiList = si.getSaiList();
    if(resAList.size() != resBList.size())
    {
        cerr << "residue list size not equal: " << resAList.size() << " " << resBList.size() << endl;
        exit(1);
    }

    float chi2A, chi2B, chi2Diff;

    for(size_t i=0;i<resAList.size();i++)
    {
        Residue* resA = resAList.at(i);
        Residue* resB = resBList.at(i);
        if(resA->altLoc != ' ') continue;
        if(resB->altLoc != ' ') continue;
        if(resA->triName != resB->triName)
        {
            cerr << "Residue type not consistent: " << resA->resID << " " << resA->triName << " " << resB->triName << endl;
            continue;
        }

        if(!resA->sidechainComplete(atLib) || !resB->sidechainComplete(atLib))
            continue;
        if(!resA->hasThreeCoreAtoms() || !resB->hasThreeCoreAtoms())
            continue;
        if(!resA->hasAtom("O")) continue;
        if(!resB->hasAtom("O")) continue;

        int chiNum = rn.chiNum(resA->triName);
        if(chiNum < 1)
            continue;
        vector<string>* scAtomNames = atLib->getAminoAcidSidechainAtomNames(resA->triName);
        chi2A = dihedral(resA->getAtom("CA")->getCoord(), resA->getAtom("CB")->getCoord(), resA->getAtom(scAtomNames->at(1))->getCoord(), resA->getAtom(scAtomNames->at(2))->getCoord());
        chi2B = dihedral(resB->getAtom("CA")->getCoord(), resB->getAtom("CB")->getCoord(), resB->getAtom(scAtomNames->at(1))->getCoord(), resB->getAtom(scAtomNames->at(2))->getCoord());
        chi2Diff = chi2A - chi2B;
        if(chi2Diff > 180)
            chi2Diff -= 360;
        if(chi2Diff < -180)
            chi2Diff += 360;
        chi2DiffList.push_back(abs(chi2Diff));
        saiList.push_back(pdbSaiList.at(i));
        resTypeList.push_back(resA->triName);
    }
}

void disulfideBondDiff(PDB* pdbA, PDB* pdbB, vector<string>& inBoth, vector<string>& onlyA, vector<string>& onlyB)
{

    vector<int> cysIndexList;
    vector<Residue*> tmpAList = pdbA->getValidResList();
    vector<Residue*> tmpBList = pdbB->getValidResList();
    if(tmpAList.size() != tmpBList.size())
    {
        cout << "list size not equal: " << pdbA->getPDBID() << endl;
        return ;
    }

    int num = 0;
    vector<Residue*> resAList;
    vector<Residue*> resBList;

    for(size_t i=0;i<tmpAList.size();i++)
    {
        if(tmpAList.at(i)->hasThreeCoreAtoms() )
        {
            resAList.push_back(tmpAList.at(i));
            resBList.push_back(tmpBList.at(i));
        }
    }



    int resNum = resAList.size();
    for(int i=0;i<resNum;i++)
    {

        if(resAList.at(i)->triName == "CYS" && resAList.at(i)->hasAtom("SG"))
            cysIndexList.push_back(i);
    }

    bool dsA, dsB;
    for(size_t i=0;i<cysIndexList.size();i++)
    {
        int indexI = cysIndexList.at(i);
        Atom* sg1A = resAList.at(indexI)->getAtom("SG");
        Atom* sg1B = resBList.at(indexI)->getAtom("SG");
        if(sg1A == NULL || sg1B == NULL)
            return;
        for(size_t j=i+1;j<cysIndexList.size();j++)
        {
            int indexJ = cysIndexList.at(j);
            Atom* sg2A = resAList.at(indexJ)->getAtom("SG");
            Atom* sg2B = resBList.at(indexJ)->getAtom("SG");
            if(sg2A == NULL || sg2B == NULL)
                return;
            if(sg1A->distance(*sg2A) < 2.5)
                dsA = true;
            else
                dsA = false;

            if(sg1B->distance(*sg2B) < 2.5)
                dsB = true;
            else
                dsB = false;


            string tag = resAList.at(indexI)->chainID+"-"+resAList.at(indexI)->resID + "-"+resAList.at(indexJ)->resID;

            if(dsA && dsB)
                inBoth.push_back(tag);
            else if(dsA)
                onlyA.push_back(tag);
            else if(dsB)
                onlyB.push_back(tag);
        }
    }
}

/*
void hydrogenBondDiff(PDB* pdbA, PDB* pdbB, vector<string>& inBoth, vector<string>& onlyA, vector<string>& onlyB)
{

    vector<Residue*> resAList = pdbA->getValidResList();
    vector<Residue*> resBList = pdbB->getValidResList();

    if(resAList.size() != resBList.size())
    {
        cerr << "residue list size not equal: " << resAList.size() << " " << resBList.size() << endl;
        exit(1);
    }

    vector<vector<Atom*>> bbPolarAtomList;
    vector<vector<Atom*>> scPolarAtomList;

    int resNum = resAList.size();
    for(int i=0;i<resNum;i++)
    {
        if(!resAList.at(i)->sidechainComplete(atLib) || !resBList.at(i)->sidechainComplete(atLib))
        {
            cerr << "sidechain not complete: " << resAList.at(i) << endl;
            exit(0);
        }



    }





}
*/

/*
void hydrogenBondCount(PDB* pdb, AtomLib* atLib, int* count, set<string>& hbs)
{
    ResName rn;
    StructureInfo si(pdb);
    si.updateTorsion();
    SasaPSD rsp;
    si.updateSAI(&rsp);

    vector<Residue*> resList = pdb->getValidResList();
    vector<pResidue*> pResList;
    vector<pScConformer*> pScList;

    map<string, int> donorTypeToIndex;
    map<string, int> acceptorTypeToIndex;

    for(int i=0;i<20;i++)
    {
        string triName = rn.intToTri(i);
        donorTypeToIndex[triName+"-N"] = 0;
        acceptorTypeToIndex[triName+"-O"] = 0;
    }

    donorTypeToIndex["HIS-ND1"] = 1;
    donorTypeToIndex["HIS-NE2"] = 1;
    donorTypeToIndex["ARG-NH1"] = 2;
    donorTypeToIndex["ARG-NH2"] = 2;
    donorTypeToIndex["ARG-NE"] = 3;
    donorTypeToIndex["LYS-NZ"] = 4;
    donorTypeToIndex["TRP-NE1"] = 5;
    donorTypeToIndex["THR-OG1"] = 6;
    donorTypeToIndex["SER-OG"] = 6;
    donorTypeToIndex["TYR-OH"] = 7;
    donorTypeToIndex["GLN-NE2"] = 8;
    donorTypeToIndex["ASN-ND2"] = 8;

    acceptorTypeToIndex["GLU-OE1"] = 1;
    acceptorTypeToIndex["GLU-OE2"] = 1;
    acceptorTypeToIndex["ASP-OD1"] = 1;
    acceptorTypeToIndex["ASP-OD2"] = 1;
    acceptorTypeToIndex["GLN-OE1"] = 2;
    acceptorTypeToIndex["ASN-OD1"] = 2;
    acceptorTypeToIndex["SER-OG"] = 3;
    acceptorTypeToIndex["THR-OG1"] = 3;
    acceptorTypeToIndex["TYR-OH"] = 4;
    acceptorTypeToIndex["HIS-ND1"] = 5;
    acceptorTypeToIndex["HIS-NE2"] = 5;

    float distAve[54];
    distAve[ 0] = 2.9300;
    distAve[ 1] = 2.8900;
    distAve[ 2] = 2.9000;
    distAve[ 3] = 3.0500;
    distAve[ 4] = 2.9800;
    distAve[ 5] = 3.0000;
    distAve[ 6] = 2.8100;
    distAve[ 7] = 2.7400;
    distAve[ 8] = 2.8000;
    distAve[ 9] = 2.7500;
    distAve[10] = 2.7100;
    distAve[11] = 3.0900;
    distAve[12] = 2.9000;
    distAve[13] = 2.8800;
    distAve[14] = 2.8900;
    distAve[15] = 2.9400;
    distAve[16] = 2.9500;
    distAve[17] = 3.0400;
    distAve[18] = 2.8600;
    distAve[19] = 2.8500;
    distAve[20] = 2.8400;
    distAve[21] = 2.9000;
    distAve[22] = 3.0800;
    distAve[23] = 2.9800;
    distAve[24] = 2.8200;
    distAve[25] = 2.7800;
    distAve[26] = 2.8100;
    distAve[27] = 2.8600;
    distAve[28] = 2.9400;
    distAve[29] = 2.8900;
    distAve[30] = 2.8900;
    distAve[31] = 2.8600;
    distAve[32] = 2.8800;
    distAve[33] = 2.9100;
    distAve[34] = 2.9300;
    distAve[35] = 2.9900;
    distAve[36] = 2.7300;
    distAve[37] = 2.6600;
    distAve[38] = 2.6900;
    distAve[39] = 2.7400;
    distAve[40] = 2.7100;
    distAve[41] = 2.7600;
    distAve[42] = 2.6700;
    distAve[43] = 2.6100;
    distAve[44] = 2.6600;
    distAve[45] = 2.7100;
    distAve[46] = 2.7100;
    distAve[47] = 2.7100;
    distAve[48] = 2.9300;
    distAve[49] = 2.9200;
    distAve[50] = 2.9300;
    distAve[51] = 2.9500;
    distAve[52] = 2.9800;
    distAve[53] = 3.0200;


    for(int i=0;i<resList.size();i++)
    {
        Residue* res = resList.at(i);
        if(!res->hasThreeCoreAtoms())
            continue;
        if(res->triName == "MSE")
            res->triName = "MET";
        if(!rn.isStandardAminoAcid(res->triName))
            continue;
        if(!res->sidechainComplete(atLib))
            continue;

        pResidue *pRes;

        if(i==0)
            pRes = new pResidue(resList.at(i),&si, atLib);
        else
            pRes = new pResidue(resList.at(i), resList.at(i-1), &si, atLib);


        Rotamer* rot = resList.at(i)->natRotamer(atLib);
        pScConformer* pSc = new pScConformer(pRes,rot,atLib);
        pResList.push_back(pRes);
        pScList.push_back(pSc);
    }


    int resN = pResList.size();

//    cout << "resNum: " << resN;
    for(int i=0;i<resN;i++)
    {
        pResidue* pResI = pResList.at(i);
        pScConformer* pScI = pScList.at(i);

        for(int j=i+2;j<resN;j++)
        {
            pResidue* pResJ = pResList.at(j);
            pScConformer* pScJ = pScList.at(j);

            if(pResI->getCB().distance(pResJ->getCB()) > 12)
                continue;

            //cout << "scI to bbJ" << endl;

            char resTag[40];
            sprintf(resTag,"%d-%d-",i,j);
            string tag(resTag);

            //cout << tag << endl;

            for(int x = 0;x<pScI->polarAtomNum;x++)
            {
                PolarAtom* paI = pScI->polarAtoms.at(x);
                for(int y=0;y<pResJ->polarAtomNum;y++)
                {
                    PolarAtom* paJ = pResJ->polarAtoms.at(y);
                    if(paI->isDonerAtom() && paJ->isAcceptorAtom())
                    {
                        int donerID = donorTypeToIndex[paI->getName()];
                        int acceptID = acceptorTypeToIndex[paJ->getName()];
                        int pairType = donerID*6+acceptID;
                        float d0Av = distAve[pairType];
                        if(paI->getCore().distance(paJ->getCore()) < d0Av + 0.3)
                        {
                            string atomTag = paI->getName()+"-"+paJ->getName();
                            //cout << "atomTag " << atomTag << endl;
                            hbs.insert(tag+"scbb");
                            count[pairType]++;
                        }
                    }
                    else if(paI->isAcceptorAtom() && paJ->isDonerAtom())
                    {
                        int acceptID = acceptorTypeToIndex[paI->getName()];
                        int donerID = donorTypeToIndex[paJ->getName()];
                        int pairType = donerID*6+acceptID;
                        float d0Av = distAve[pairType];
                        if(paI->getCore().distance(paJ->getCore()) < d0Av + 0.3)
                        {
                            string atomTag = paI->getName()+"-"+paJ->getName();
                            hbs.insert(tag+"scbb");
                            count[pairType]++;
                        }
                    }


                }
            }

            //cout << "bbI to scJ" << endl;
            for(int x = 0;x<pResI->polarAtomNum;x++)
            {
                PolarAtom* paI = pResI->polarAtoms.at(x);
                for(int y=0;y<pScJ->polarAtomNum;y++)
                {
                    PolarAtom* paJ = pScJ->polarAtoms.at(y);
                    if(paI->isDonerAtom() && paJ->isAcceptorAtom())
                    {
                        int donerID = donorTypeToIndex[paI->getName()];
                        int acceptID = acceptorTypeToIndex[paJ->getName()];
                        int pairType = donerID*6+acceptID;
                        float d0Av = distAve[pairType];
                        if(paI->getCore().distance(paJ->getCore()) < d0Av + 0.3)
                        {
                            string atomTag = paI->getName()+"-"+paJ->getName();
                            hbs.insert(tag+"bbsc");
                            count[pairType]++;
                        }
                    }
                    else if(paI->isAcceptorAtom() && paJ->isDonerAtom())
                    {
                        int acceptID = acceptorTypeToIndex[paI->getName()];
                        int donerID = donorTypeToIndex[paJ->getName()];
                        int pairType = donerID*6+acceptID;
                        float d0Av = distAve[pairType];
                        if(paI->getCore().distance(paJ->getCore()) < d0Av + 0.3)
                        {
                            string atomTag = paI->getName()+"-"+paJ->getName();
                            hbs.insert(tag+"bbsc");
                            count[pairType]++;
                        }
                    }

                }
            }



        //    cout << "scI to scJ" << endl;
            for(int x = 0;x<pScI->polarAtomNum;x++)
            {
                PolarAtom* paI = pScI->polarAtoms.at(x);
                for(int y=0;y<pScJ->polarAtomNum;y++)
                {
                    PolarAtom* paJ = pScJ->polarAtoms.at(y);
                    if(paI->isDonerAtom() && paJ->isAcceptorAtom())
                    {
                        int donerID = donorTypeToIndex[paI->getName()];
                        int acceptID = acceptorTypeToIndex[paJ->getName()];
                        int pairType = donerID*6+acceptID;
                        float d0Av = distAve[pairType];
                        if(paI->getCore().distance(paJ->getCore()) < d0Av + 0.3)
                        {
                            string atomTag = paI->getName()+"-"+paJ->getName();
                            hbs.insert(tag+"scsc");
                            count[pairType]++;
                        }
                    }
                    else if(paI->isAcceptorAtom() && paJ->isDonerAtom())
                    {
                        int acceptID = acceptorTypeToIndex[paI->getName()];
                        int donerID = donorTypeToIndex[paJ->getName()];
                        int pairType = donerID*6+acceptID;
                        float d0Av = distAve[pairType];
                        if(paI->getCore().distance(paJ->getCore()) < d0Av + 0.3)
                        {
                            string atomTag = paI->getName()+"-"+paJ->getName();
                            hbs.insert(tag+"scsc");
                            count[pairType]++;
                        }
                    }

                }
            }

        //    cout << "finish" << endl;
        }
    }
}
*/


} /* namespace NSPdesignseq */

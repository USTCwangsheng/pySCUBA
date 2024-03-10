/*
 * SeqMinimize.h
 *
 *  Created on: 2017年12月26日
 *      Author: notxp
 */

#ifndef DESIGNSEQ_SEQMINIMIZE_H_
#define DESIGNSEQ_SEQMINIMIZE_H_

#include <string>
#include <vector>
#include "stdio.h"
#include <iostream>
#include "time.h"
#include "designseq/RotamerGroup.h"
#include "designseq/DesignTemplate.h"
#include "designseq/PackingTemplate.h"
#include "designseq/SingleSitePackingTemplate.h"

namespace NSPdesignseq {

using namespace std;



/*
 * added by xuyang, 2019.4.4
 */
struct AANumberControlInMC {
    std::set<int> sites;
    char aa;
    int min;
    int max;
    int num;
};





class RotSequence {
private:
    int seqLength;
    int* seqChoices;
    char* aaSeq;
    vector<RotamerGroup*> rotGroups;
    ResName rn;

public:
    RotSequence(int len);
    RotSequence(vector<RotamerGroup*>& groups);
    void copyValue(RotSequence* from);
    void setRandomChoice();
    void applyMutation(int pos, int choice);

    string toAASequence();
    int getChoice(int pos) { return this->seqChoices[pos];}
    Rotamer* getRotamer(int pos, int choice) {return this->rotGroups[pos]->rotList[choice];}
    Rotamer* getRotamer(int pos) {
        int choice = this->seqChoices[pos];
        return this->rotGroups[pos]->rotList[choice];
    }
    int getLength() {return this->seqLength;}
    int choiceNum(int pos) {return this->rotGroups[pos]->rotNum;}

    float maxIdentityTo(vector<string>& preSeqs);
    float mutIdentityEnergy(vector<string>& preSeqs, int pos, int mutChoice, double idtCutoff, double sigma);
    float mutNativeSeqIdtEnergy(const string& natSeq, int pos, int mutChoice, double idtCutoff, double sigma);

    float identity(RotSequence* other) {
        int choiceA, choiceB;
        string typeA, typeB;
        int idtCount = 0;
        for(int pos=0;pos<seqLength;pos++){
            choiceA = this->seqChoices[pos];
            choiceB = other->seqChoices[pos];
            if(choiceA == choiceB){
                idtCount++;
                continue;
            }
            typeA = this->rotGroups[pos]->rotList[choiceA]->triName;
            typeB = other->rotGroups[pos]->rotList[choiceB]->triName;
            if(typeA == typeB)
                idtCount++;
        }
        return 1.0*idtCount/seqLength;
    }

    int idtCount(RotSequence* other) {
        int choiceA, choiceB;
        string typeA, typeB;
        int idtCount = 0;
        for(int pos=0;pos<seqLength;pos++){
            choiceA = this->seqChoices[pos];
            choiceB = other->seqChoices[pos];
            if(choiceA == choiceB){
                idtCount++;
                continue;
            }
            typeA = this->rotGroups[pos]->rotList[choiceA]->triName;
            typeB = other->rotGroups[pos]->rotList[choiceB]->triName;
            if(typeA == typeB)
                idtCount++;
        }
        return idtCount;
    }

    RotamerGroup* getRotGroup(int pos) {return this->rotGroups[pos];}

    double mutEnergy(int pos, int mutChoice, EnergyCalculatorTemplate* ec);
    double resInvolvedEnergy(int pos, int choice, EnergyCalculatorTemplate* ec);
    double resInvolvedHbondEnergy(int pos, SingleSitePackingTemplate* st);
    double pairInvolvedEnergy(int posA, int choiceA, int posB, int choiceB, EnergyCalculatorTemplate* ec);
    void acceptMut(int pos, int mutChoice);
    void acceptMutAndUpdateSeq(int pos, int mutChoice);
    double totEnergy(EnergyCalculatorTemplate* ec);
    void calcEnergyComponents(const DesignTemplate* dt, double* s1, double* rot, double* ref) const;
    void printDetailEnergy(EnergyCalculatorTemplate* ec);

    void checkChoice(){
        for(int i=0;i<seqLength;i++){
            RotamerGroup* group = rotGroups.at(i);
            int rotNum = group->rotNum;
            int choice = this->seqChoices[i];
            cout << "site: " << i << " " << rotNum << " " << choice << endl;
        }
    }



    /*
     * added by xuyang, 2019.3.31
     */

    char getaatype(int npos) {return aaSeq[npos];}
    ResName & getrn() {return rn;}
    RotSequence(vector<RotamerGroup*>& groups, std::vector<AANumberControlInMC> &aacmc);
    bool adjustSeq_1time(AANumberControlInMC &ac);
    std::vector<int> goodChoiceOnly(char aa, int pos);
    std::vector<int> goodChoiceExcept(char aa, int pos);


    virtual ~RotSequence();
};



class DesignMC {
private:
    EnergyCalculatorTemplate* dt;
    double T0;
    double T1;
    double annealFactor;
    int step;

public:
    DesignMC(DesignTemplate* dt){
        this->dt = dt;
        this->T0 = 10.0;
        this->T1 = 0.001;
        this->annealFactor = 0.95;
        this->step = 50000;
        srand((unsigned)time(NULL));
    }

    DesignMC(PackingTemplate* pt){
        this->dt = pt;
        this->T0 = 10.0;
        this->T1 = 0.01;
        this->annealFactor = 0.9;
        this->step = 10000;
        srand((unsigned)time(NULL));
    }

    DesignMC(SingleSitePackingTemplate* st){
        this->dt = st;
        this->T0 = 5.0;
        this->T1 = 0.001;
        this->annealFactor = 0.9;
        this->step = 10000;
        srand((unsigned)time(NULL));
    }

    bool accept(double mutEnergy, double T);


    int minEnergyChoice(RotSequence* unit, int pos);
    pair<int,int> minEnergyPairChoice(RotSequence* unit, int posA, int posB);

    double mcRun(RotSequence* result);
    double mcRunWithIdtRestrainToNat(RotSequence* result, const string& natSeq, double idtCutoff, double sigma);
    double mcRunWithIdtRestrain(RotSequence* result, vector<string>& preSeqs, double idtCutoff, double sigma);

    void packing(RotSequence* result);
    void packing2(RotSequence* result);
    void packing3(RotSequence* result);

    void printPDB(RotSequence* unit, string outputFile);



    /*
     * added by xuyang, 2019.3.31
     */
    double mcRunWithAACountRestrain(RotSequence* result, std::vector<AANumberControlInMC> &aacmc);
    double mcRunWithIdtAACountRestrain(RotSequence* result, vector<string>& preSeqs,
                double idtCutoff,double sigma, std::vector<AANumberControlInMC> &aacmc);


    virtual ~DesignMC();
};




class DoubleTargetDesign{
private:
    DesignTemplate* dt1;
    DesignTemplate* dt2;

    vector<RotamerGroup*> rotGroups;
    vector<vector<int>> posAAList;
    vector<vector<vector<int>>> posAAChoiceList;


    double weight1;
    double weight2;

    double T0;
    double T1;
    double annealFactor;
    int step;
public:
    DoubleTargetDesign(DesignTemplate* dt1, DesignTemplate* dt2, double wt1, double wt2);

    void loadEnergyTable(S1EnergyTable& s1ET, S2MatrixFinder& s2ET);
    void mcRun(RotSequence* result1, RotSequence* result2);

    void doubleTargetDesignMCIdt(RotSequence* result1, RotSequence* result2, double idt, double idtFactor);
    void threeTargetDesignMC(RotSequence* result1, RotSequence* result2, RotSequence* result3, RotSequence* result4, double idt, double idtFactor);
    void threeTargetDesignMC34Fixed(RotSequence* result1, RotSequence* result2, RotSequence* result3, RotSequence* result4, double idt, double idtFactor);


    bool accept(double mutEnergy, double T);
    void printPDB(RotSequence* unit, string outputFile, DesignTemplate* dt);
    virtual ~DoubleTargetDesign();
};



} /* namespace NSPdesignseq */

#endif /* DESIGNSEQ_SEQMINIMIZE_H_ */

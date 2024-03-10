/*
 * SeqMinimize.cpp
 *
 *  Created on: 2017��12��26��
 *      Author: notxp
 */

#include "designseq/SeqMinimize.h"

namespace NSPdesignseq {



    RotSequence::RotSequence(int len){
        this->seqLength = len;
        this->seqChoices = new int[len];
        this->aaSeq = new char[len];
        for(int i=0;i<len;i++) {
            this->seqChoices[i] = 0;
            this->aaSeq[i] = 'A';
        }

    }

    RotSequence::RotSequence(vector<RotamerGroup*>& groups){

        int len = groups.size();
        this->seqLength = len;
        this->seqChoices = new int[len];
        this->aaSeq = new char[len];
        for(int i=0;i<len;i++){
            this->rotGroups.push_back(groups[i]);
        }
        setRandomChoice();
    }

    void RotSequence::copyValue(RotSequence* from){
        if(this->seqLength != from->seqLength){
            cerr << "rotamer sequence length not equal" << endl;
            exit(1);
        }

        this->rotGroups.clear();
        for(int i=0;i<seqLength;i++){
            this->seqChoices[i] = from->seqChoices[i];
            this->aaSeq[i] = from->aaSeq[i];
            this->rotGroups.push_back(from->rotGroups[i]);
        }
    }

    void RotSequence::setRandomChoice(){
    //    cout << "set Random choice:" << endl;
        for(int i=0;i<this->seqLength;i++){
            int choiceNum = this->rotGroups[i]->rotNum;
    //        cout << "rot number: " << choiceNum << endl;
            int randNum = rand()%choiceNum;
            this->seqChoices[i] = randNum;
            Rotamer* rot = this->rotGroups[i]->rotList[this->seqChoices[i]];
            aaSeq[i] = rn.triToSin(rot->triName);
        }
    }

    void RotSequence::applyMutation(int pos, int choice){
        this->seqChoices[pos] = choice;
        this->aaSeq[pos] = rn.triToSin(this->rotGroups[pos]->rotList[choice]->triName);
    }

    string RotSequence::toAASequence(){
        ResName rn;
        char seq[this->seqLength+1];
        for(int i=0;i<seqLength;i++){
            Rotamer* rot = this->rotGroups[i]->rotList[this->seqChoices[i]];
            seq[i] = rn.triToSin(rot->triName);
        }
        seq[this->seqLength] = '\0';
        string s = seq;
        return s;
    }

    float RotSequence::maxIdentityTo(vector<string>& preSeqs){
        int maxIdtCount = 0;
        int n = preSeqs.size();
        for(int i=0;i<n;i++){
            int count = 0;
            for(int j=0;j<seqLength;j++){
                if(aaSeq[j] == preSeqs[i][j])
                    count ++;
            }
            if(count > maxIdtCount){
                maxIdtCount = count;
            }
        }
        return maxIdtCount*1.0/seqLength;
    }

    float RotSequence::mutIdentityEnergy(vector<string>& preSeqs, int pos, int mutChoice, double idtCutoff, double sigma){
        int bMax = 0;
        int aMax = 0;
        int n = preSeqs.size();
        char mutAAType = rn.triToSin(this->rotGroups[pos]->rotList[mutChoice]->triName);
        for(int i=0;i<n;i++){
            int beforeMutCount = 0;
            int afterMutCount = 0;
            for(int j=0;j<seqLength;j++){
                if(j == pos) {
                    if(mutAAType == preSeqs[i][j])
                        afterMutCount ++;
                    if(aaSeq[j] == preSeqs[i][j])
                        beforeMutCount++;
                }
                else if(aaSeq[j] == preSeqs[i][j]) {
                    beforeMutCount++;
                    afterMutCount ++;
                }
            }
            if(beforeMutCount > bMax){
                bMax = beforeMutCount;
            }
            if(afterMutCount > aMax) {
                aMax = afterMutCount;
            }
        }

        double beforeMutIdt = 1.0*bMax/seqLength;
        double afterMutIdt = 1.0*aMax/seqLength;
        double e1 = exp((beforeMutIdt-idtCutoff)/sigma);
        double e2 = exp((afterMutIdt-idtCutoff)/sigma);
        return e2-e1;
    }

    float RotSequence::mutNativeSeqIdtEnergy(const string& natSeq, int pos, int mutChoice, double idtCutoff, double sigma){
        char mutAAType = rn.triToSin(this->rotGroups[pos]->rotList[mutChoice]->triName);
        int beforeMutCount = 0;
        int afterMutCount = 0;
        for(int j=0;j<seqLength;j++){
            if(j == pos) {
                if(mutAAType == natSeq[j])
                    afterMutCount ++;
                if(aaSeq[j] == natSeq[j])
                    beforeMutCount++;
            }
            else if(aaSeq[j] == natSeq[j]) {
                beforeMutCount++;
                afterMutCount ++;
            }
        }
        double beforeMutIdt = 1.0*beforeMutCount/seqLength;
        double afterMutIdt = 1.0*afterMutCount/seqLength;
        double e1 = exp((beforeMutIdt-idtCutoff)*(beforeMutIdt-idtCutoff)/sigma);
        double e2 = exp((afterMutIdt-idtCutoff)*(afterMutIdt-idtCutoff)/sigma);
        return e2-e1;

    }


    double RotSequence::mutEnergy(int pos, int mutChoice, EnergyCalculatorTemplate* ec){
        int oldChoice = this->seqChoices[pos];
        if(mutChoice == oldChoice) return 0;
        return resInvolvedEnergy(pos, mutChoice, ec) - resInvolvedEnergy(pos, oldChoice, ec);
    }



    double RotSequence::resInvolvedEnergy(int pos, int choice, EnergyCalculatorTemplate* ec){
        double e1 = ec->eaList[pos]->getEnergy(choice);
        double e2 = 0;
        int pairNum = ec->involvedPairs[pos]->size();
        for(int i=0;i<pairNum;i++){
            int pairID = ec->involvedPairs[pos]->at(i);
            BackboneSitesPair* pair = ec->bsPairs[pairID];
            EnergyMatrix* em = ec->emList[pairID];
            int posA = pair->siteA->resseq;
            int posB = pair->siteB->resseq;
            if(pos == posA){
                e2 += em->getEnergy(choice,this->seqChoices[posB]);
            }
            else if(pos == posB){
                e2 += em->getEnergy(this->seqChoices[posA], choice);
            }
            else{
                cerr << "posA: " << posA << " posB: " << posB << " mutPos: " << pos << " resInvolved Energy error" << endl;
                exit(1);
            }
        }
        return e1+e2;
    }

    double RotSequence::resInvolvedHbondEnergy(int pos, SingleSitePackingTemplate* st) {
        int pairNum = st->involvedPairs[pos]->size();
        double hbEnergy = 0;
        for(int i=0;i<pairNum;i++){
            int pairID = st->involvedPairs[pos]->at(i);
            BackboneSitesPair* pair = st->bsPairs[pairID];
            int posA = pair->siteA->resseq;
            int posB = pair->siteB->resseq;
            if(abs(posB-posA) == 1) continue;
            int choiceA = this->seqChoices[posA];
            int choiceB = this->seqChoices[posB];
            Rotamer* rotA = st->rotGroups[posA]->rotList[choiceA];
            Rotamer* rotB = st->rotGroups[posB]->rotList[choiceB];
            Conformer* confA = new Conformer(rotA, st->bsList[posA], &st->atLib);
            Conformer* confB = new Conformer(rotB, st->bsList[posB], &st->atLib);
            int index = 0;
            if(pos == posB) index = 1;
            hbEnergy += pairHBEnergyPack(confA, confB, index, &st->aec);
            delete confA;
            delete confB;
        }
        return hbEnergy;
    }

    double RotSequence::pairInvolvedEnergy(int posX, int choiceX, int posY, int choiceY, EnergyCalculatorTemplate* ec){
        double e1 = ec->eaList[posX]->getEnergy(choiceX) + ec->eaList[posY]->getEnergy(choiceY);
        double e2 = 0;
        int pairNumA = ec->involvedPairs[posX]->size();
        for(int i=0;i<pairNumA;i++){
            int pairID = ec->involvedPairs[posX]->at(i);
            BackboneSitesPair* pair = ec->bsPairs[pairID];
            EnergyMatrix* em = ec->emList[pairID];
            int posA = pair->siteA->resseq;
            int posB = pair->siteB->resseq;
            if(posX == posA && posY == posB)
                e2 += em->getEnergy(choiceX,choiceY);
            else if(posX == posB && posY == posA)
                e2 += em->getEnergy(choiceY,choiceX);
            else if(posX == posA){
                e2 += em->getEnergy(choiceX,this->seqChoices[posB]);
            }
            else if(posX == posB){
                e2 += em->getEnergy(this->seqChoices[posA], choiceX);
            }
            else{
                cerr << "posA: " << posA << " posB: " << posB << " mutPos: " << posX << " resInvolved Energy error" << endl;
                exit(1);
            }
        }

        int pairNumB = ec->involvedPairs[posY]->size();
        for(int i=0;i<pairNumB;i++){
            int pairID = ec->involvedPairs[posY]->at(i);
            BackboneSitesPair* pair = ec->bsPairs[pairID];
            EnergyMatrix* em = ec->emList[pairID];
            int posA = pair->siteA->resseq;
            int posB = pair->siteB->resseq;
            if((posX == posA && posY == posB) || (posX == posB && posY == posA))
                continue;
            if(posY == posA){
                e2 += em->getEnergy(choiceY,this->seqChoices[posB]);
            }
            else if(posY == posB){
                e2 += em->getEnergy(this->seqChoices[posA], choiceY);
            }
            else{
                cerr << "posA: " << posA << " posB: " << posB << " mutPos: " << posY << " resInvolved Energy error" << endl;
                exit(1);
            }
        }
        return e1+e2;
    }

    void RotSequence::acceptMut(int pos, int mutChoice){
        this->seqChoices[pos] = mutChoice;
    }

    void RotSequence::acceptMutAndUpdateSeq(int pos, int mutChoice){
        this->seqChoices[pos] = mutChoice;
        this->aaSeq[pos] = rn.triToSin(this->rotGroups[pos]->rotList[mutChoice]->triName);
    }

    double RotSequence::totEnergy(EnergyCalculatorTemplate* ec){
        double e = 0;
        for(int i=0;i<ec->eaList.size();i++){

            e += ec->eaList.at(i)->getEnergy(this->seqChoices[i]);
        }
        for(int i=0;i<ec->emList.size();i++){
            BackboneSitesPair* pair = ec->bsPairs[i];
            EnergyMatrix* em = ec->emList[i];
            int posA = pair->siteA->resseq;
            int posB = pair->siteB->resseq;
            e += em->getEnergy(this->seqChoices[posA], this->seqChoices[posB]);
        }
        return e;
    }

    void RotSequence::printDetailEnergy(EnergyCalculatorTemplate* ec){
        double e = 0;
        double e1 = 0;
        double e2 = 0;
        for(int i=0;i<ec->eaList.size();i++){
            Rotamer* rot = getRotamer(i);
            e = ec->eaList[i]->getEnergy(getChoice(i));
            e1 += e;
            printf("Res: %3d Rot: %7s Energy: %10.3f\n", i+1, rot->rotName.c_str(), e);
        }

        int n = ec->emList.size();
        for(int i=0;i<n;i++){
            BackboneSitesPair* pair = ec->bsPairs[i];
            EnergyMatrix* em = ec->emList[i];
            int posA = pair->siteA->resseq;
            int posB = pair->siteB->resseq;
            Rotamer* rotA = getRotamer(posA);
            Rotamer* rotB = getRotamer(posB);
            e = em->getEnergy(this->seqChoices[posA], this->seqChoices[posB]);
            e2+= e;
            printf("Pair: %3d %3d Rot: %-7s %-7s Energy: %10.3f\n", posA+1, posB+1, rotA->rotName.c_str(), rotB->rotName.c_str(), e );
        }

        printf("TotalEnergy: E1: %7.3f E2: %7.3f\n", e1, e2);
    }

    RotSequence::~RotSequence(){
        delete[] this->seqChoices;
    }

    bool DesignMC::accept(double mutEnergy, double T){
        if(T == 9){
            if(mutEnergy > 0)
                return false;
            else
                return true;
        }
        else if(mutEnergy > 0){
            float pAc = exp(-1.0*mutEnergy/T);
            float r = 1.0*rand()/RAND_MAX;
            return r < pAc;
        }
        else
            return true;
    }

    void DesignMC::packing(RotSequence* result){


        RotSequence unit(dt->rotGroups);
        result->copyValue(&unit);
        int len = unit.getLength();

        int randPos, randChoice;
        double mutEnergy;

        double energy = 0;
        double minEnergy = energy;
        int count = 0;
        step = len*300;
        for(double T = 5.0;T>0.01;T=T*annealFactor){

            bool noChange = true;
            int acNum = 0;
            for(int i=0;i<step;i++){
                count ++;
                randPos = rand()%len;
                randChoice = rand()%(unit.choiceNum(randPos));
                mutEnergy = unit.mutEnergy(randPos,randChoice, dt);
                if(mutEnergy < 0)
                    noChange = false;
                if(accept(mutEnergy, T)){
                    acNum++;
                    unit.acceptMut(randPos, randChoice);
                    energy += mutEnergy;
                    if(energy < minEnergy){
                        minEnergy = energy;
                        result->copyValue(&unit);
                    }
                }
            }
    //        printf("%7.4f %10.4f %10.4f %4d\n",T,energy,minEnergy,acNum);
            if(noChange) break;
        }


        for(int pos=0;pos<len;pos++){
            for(int i=0;i<unit.choiceNum(pos);i++){
                mutEnergy = unit.mutEnergy(pos,i,dt);
                if(mutEnergy < 0){
                       unit.acceptMut(randPos, randChoice);
                       energy += mutEnergy;
                       if(energy < minEnergy){
                           minEnergy = energy;
                           result->copyValue(&unit);
                       }
                }
            }
        }

        for(int pos=0;pos<len;pos++){
            for(int i=0;i<unit.choiceNum(pos);i++){
                mutEnergy = unit.mutEnergy(pos,i,dt);
                if(mutEnergy < 0){
                       unit.acceptMut(randPos, randChoice);
                       energy += mutEnergy;
                       if(energy < minEnergy){
                           minEnergy = energy;
                           result->copyValue(&unit);
                       }
                }
            }
        }

    }

    void DesignMC::packing2(RotSequence* result){
        RotSequence unit(dt->rotGroups);
        result->copyValue(&unit);
        int len = unit.getLength();
        int randPos, randChoice;
        double mutEnergy;
        double energy = 0;
        double minEnergy = energy;
        int count = 0;
        step = len*1000;
       for(double T = this->T0;T>this->T1;T=T*annealFactor){
         bool noChange = true;
            int acNum = 0;
            for(int i=0;i<step;i++){
                count ++;
                randPos = rand()%len;
                randChoice = rand()%(unit.choiceNum(randPos));
                mutEnergy = unit.mutEnergy(randPos,randChoice, dt);
                if(mutEnergy < 0)
                    noChange = false;
                if(accept(mutEnergy, T)){
                    acNum++;
                    unit.acceptMut(randPos, randChoice);
                    energy += mutEnergy;
                    if(energy < minEnergy){
                        minEnergy = energy;
                        result->copyValue(&unit);
                    }
                }
            }
         }

        int pairNum = dt->bsPairs.size();
        for(int i=0;i<pairNum;i++){
            int posA = dt->bsPairs[i]->siteA->resseq;
            int posB = dt->bsPairs[i]->siteB->resseq;
            pair<int,int> p = minEnergyPairChoice(result,posA,posB);
            result->acceptMut(posA,p.first);
            result->acceptMut(posB,p.second);
        }

        for(int i=0;i<len;i++){
            int c = minEnergyChoice(result, i);
            result->acceptMut(i,c);
        }

        for(int i=0;i<pairNum;i++){
            int posA = dt->bsPairs[i]->siteA->resseq;
            int posB = dt->bsPairs[i]->siteB->resseq;
            pair<int,int> p = minEnergyPairChoice(result,posA,posB);
            result->acceptMut(posA,p.first);
            result->acceptMut(posB,p.second);
        }

        for(int i=0;i<len;i++){
            int c = minEnergyChoice(result, i);
            result->acceptMut(i,c);
        }

    }

    void DesignMC::packing3(RotSequence* result){
        RotSequence unit(dt->rotGroups);
        result->copyValue(&unit);
        int len = unit.getLength();
        int randPos, randChoice;
        double mutEnergy;
        double energy = 0;
        double minEnergy = energy;
        int count = 0;
        vector<int> mutablePosList;
        for(int i=0;i<len;i++){
            if(unit.choiceNum(i) > 1)
                mutablePosList.push_back(i);
        }
        int mutableSiteNum = mutablePosList.size();

        step = mutableSiteNum*1000;
        for(double T = 5.0;T>0.01;T=T*annealFactor){
            int acNum = 0;
            for(int i=0;i<step;i++){
                count ++;
                randPos = mutablePosList[rand()%mutableSiteNum];
                randChoice = rand()%(unit.choiceNum(randPos));
                mutEnergy = unit.mutEnergy(randPos,randChoice, dt);
                if(accept(mutEnergy, T)){
                    acNum++;
                    unit.acceptMut(randPos, randChoice);
                    energy += mutEnergy;
                    if(energy < minEnergy){
                        minEnergy = energy;
                        result->copyValue(&unit);
                    }
                }
            }
         }

        int pairNum = dt->bsPairs.size();
        for(int i=0;i<pairNum;i++){
            int posA = dt->bsPairs[i]->siteA->resseq;
            int posB = dt->bsPairs[i]->siteB->resseq;
            pair<int,int> p = minEnergyPairChoice(result,posA,posB);
            result->acceptMut(posA,p.first);
            result->acceptMut(posB,p.second);
        }

        for(int i=0;i<len;i++){
            int c = minEnergyChoice(result, i);
            result->acceptMut(i,c);
        }

        for(int i=0;i<len;i++){
            int c = minEnergyChoice(result, i);
            result->acceptMut(i,c);
        }
    }

    int DesignMC::minEnergyChoice(RotSequence* unit, int pos){
        int choiceNum = unit->choiceNum(pos);
        if(choiceNum == 1) return 0;
        float minE = 99999999999.9;
        int minChoice = 0;

        float e=0;
        for(int i=0;i<choiceNum;i++){
            e = unit->resInvolvedEnergy(pos,i,this->dt);
            if(e < minE){
                minE = e;
                minChoice = i;
            }
        }
        return minChoice;
    }

    pair<int,int> DesignMC::minEnergyPairChoice(RotSequence* unit, int posA, int posB){
        int choiceNumA = unit->choiceNum(posA);
        int choiceNumB = unit->choiceNum(posB);
        float minE = 99999999999.9;
        int minChoiceA = 0;
        int minChoiceB = 0;

        float e;
        for(int i=0;i<choiceNumA;i++){
            for(int j=0;j<choiceNumB;j++){
                e = unit->pairInvolvedEnergy(posA,i,posB,j,dt);
                if(e < minE){
                    minE = e;
                    minChoiceA = i;
                    minChoiceB = j;
                }
            }
        }
        pair<int,int> p(minChoiceA,minChoiceB);
        return p;
    }

    double DesignMC::mcRunWithIdtRestrain(RotSequence* result, vector<string>& preSeqs, double idtCutoff,double sigma){

         RotSequence unit(dt->rotGroups);
         result->copyValue(&unit);
         double energy = unit.totEnergy(dt);

         int len = unit.getLength();

         int randPos, randChoice;
         double mutEnergy;
         double idtMutEnergy;

         srand((unsigned)time(NULL));
         double minEnergy = energy;
         int count = 0;
         for(double T = T0;T>T1;T=T*annealFactor){
             bool noChange = true;
             int acNum = 0;
             for(int i=0;i<step;i++){
                 count ++;

                 randPos = rand()%len;
                 randChoice = rand()%(unit.choiceNum(randPos));
                 mutEnergy = unit.mutEnergy(randPos,randChoice, dt);
                 idtMutEnergy = unit.mutIdentityEnergy(preSeqs,randPos,randChoice,idtCutoff,sigma);
                 mutEnergy += idtMutEnergy;

                 if(mutEnergy < 0)
                     noChange = false;
                 if(accept(mutEnergy, T)){
                     acNum++;
                     unit.acceptMutAndUpdateSeq(randPos, randChoice);
                     energy += mutEnergy;
                     if(energy < minEnergy){
                         minEnergy = energy;
                         result->copyValue(&unit);
                     }
                 }

             }
             printf("%7.4f %10.4f %10.4f %4d\n",T,energy,minEnergy,acNum);
             if(noChange) break;
         }

         return minEnergy;
    }



    double DesignMC::mcRun(RotSequence* result){
         RotSequence unit(dt->rotGroups);
         result->copyValue(&unit);
         double energy = unit.totEnergy(dt);
         int len = unit.getLength();
         int randPos, randChoice;
         double mutEnergy;

         srand((unsigned)time(NULL));
         double minEnergy = energy;
         int count = 0;
         for(double T = T0;T>T1;T=T*annealFactor){
             bool noChange = true;
             int acNum = 0;
             for(int i=0;i<step;i++){
                 count ++;
                 randPos = rand()%len;
                 randChoice = rand()%(unit.choiceNum(randPos));
                 mutEnergy = unit.mutEnergy(randPos,randChoice, dt);
                 if(mutEnergy < 0)
                     noChange = false;
                 if(accept(mutEnergy, T)){
                     acNum++;
                     unit.acceptMut(randPos, randChoice);
                     energy += mutEnergy;
                     if(energy < minEnergy){
                         minEnergy = energy;
                         result->copyValue(&unit);
                     }
                 }
             }
             printf("%7.4f %10.4f %10.4f %4d\n",T,energy,minEnergy,acNum);
             if(noChange) break;
         }
         double finalEnergy = result->totEnergy(dt);

         return finalEnergy;
    }

    double DesignMC::mcRunWithIdtRestrainToNat(RotSequence* result, const string& natSeq, double idtCutoff, double sigma){
        RotSequence unit(dt->rotGroups);
        result->copyValue(&unit);
        double energy = unit.totEnergy(dt);

        int len = unit.getLength();

        int randPos, randChoice;
        double mutEnergy;
        double idtMutEnergy;

        srand((unsigned)time(NULL));
        double minEnergy = energy;
        int count = 0;
        for(double T = T0;T>T1;T=T*annealFactor){
            bool noChange = true;
            int acNum = 0;
            for(int i=0;i<step;i++){
                count ++;

                randPos = rand()%len;
                randChoice = rand()%(unit.choiceNum(randPos));
                mutEnergy = unit.mutEnergy(randPos,randChoice, dt);
                idtMutEnergy = unit.mutNativeSeqIdtEnergy(natSeq, randPos, randChoice, idtCutoff, sigma);
                mutEnergy += idtMutEnergy;

                if(mutEnergy < 0)
                    noChange = false;
                if(accept(mutEnergy, T)){
                    acNum++;
                    unit.acceptMutAndUpdateSeq(randPos, randChoice);
                    energy += mutEnergy;
                    if(energy < minEnergy){
                        minEnergy = energy;
                        result->copyValue(&unit);
                    }
                }

            }
            printf("%7.4f %10.4f %10.4f %4d\n",T,energy,minEnergy,acNum);
            if(noChange) break;
        }

        double finalEnergy = result->totEnergy(dt);
        return finalEnergy;
    }

    void DesignMC::printPDB(RotSequence* unit, string outputFile){


        ProteinChain pc;
        char s[20];
        for(int i=0;i<dt->bsList.size();i++){

            BackBoneSite* bs = dt->bsList.at(i);
            Rotamer* rot = unit->getRotamer(i);
            Residue* res = new Residue(bs->strresid, bs->chainid, rot->triName);
            res->addAtom(new Atom("N", bs->ncrd()));
            res->addAtom(new Atom("CA", bs->cacrd()));
            res->addAtom(new Atom("C", bs->ccrd()));
            res->addAtom(new Atom("O", bs->ocrd()));
            res->buildRotamer(rot);
            pc.addResidue(res);
        }

        ofstream output(outputFile, ios::out);

        if(!output.is_open())
        {
            cout << "fail to open file " << outputFile << endl;
            exit(1);
        }
        pc.printPDBFormat(output, 1);


        Residue* p1;
        Atom* p2;
        for(int i=0;i<pc.getChainLength();i++){
            p1 = pc.getResList().at(i);
            for(int j=0;j<p1->getAtomList()->size();j++){
                p2 = p1->getAtomList()->at(j);
                delete p2;
            }
            delete p1;
        }
        output.close();
    }

    DesignMC::~DesignMC(){

    }

    bool DoubleTargetDesign::accept(double mutEnergy, double T){
         if(T == 9){
             if(mutEnergy > 0)
                 return false;
             else
                 return true;
         }
         else if(mutEnergy > 0){
             float pAc = exp(-1.0*mutEnergy/T);
             float r = 1.0*rand()/RAND_MAX;
             return r < pAc;
         }
         else
             return true;
     }

    DoubleTargetDesign::DoubleTargetDesign(DesignTemplate* dt1, DesignTemplate* dt2, double wt1, double wt2){
        this->dt1 = dt1;
        this->dt2 = dt2;
        this->T0 = 100.0;
        this->T1 = 0.001;
        this->annealFactor = 0.95;
        this->step = 50000;
        srand((unsigned)time(NULL));
        this->weight1 = wt1;
        this->weight2 = wt2;

        ResName rn;
        int resNum = dt1->resNum;

        vector<RotamerGroup*> groups1 = dt1->rotGroups;
        vector<RotamerGroup*> groups2 = dt2->rotGroups;

        for(int i=0;i<groups1.size();i++){
            RotamerGroup* gp = new RotamerGroup();
            groups1[i]->merge(*groups2[i],*gp);
            this->rotGroups.push_back(gp);
        }

        dt1->setRotamerGroups(this->rotGroups);
        dt1->updateResPairs();
        dt2->setRotamerGroups(this->rotGroups);
        dt2->updateResPairs();


        for(int i=0;i<resNum;i++){
            RotamerGroup* gp = rotGroups[i];
            vector<int> aaList;
            vector<vector<int>> aaChoiceList;
            for(int j=0;j<20;j++){
                vector<int> aaChoices;
                for(int k=0;k<gp->rotList.size();k++){
                    Rotamer* rot = gp->rotList[k];
                    int aaType = rn.triToInt(rot->triName);
                    if(aaType == j)
                        aaChoices.push_back(k);
                }
                aaChoiceList.push_back(aaChoices);
                if(gp->aaRots.at(j)->size() > 0){
                    aaList.push_back(j);
                }
            }
            posAAList.push_back(aaList);
            posAAChoiceList.push_back(aaChoiceList);
        }
    }

    void DoubleTargetDesign::loadEnergyTable(S1EnergyTable& s1ET, S2MatrixFinder& s2ET){
        this->dt1->loadS1S2(s1ET, s2ET);
        this->dt2->loadS1S2(s1ET, s2ET);

        this->dt1->loadSingleResidueEnergy();
        this->dt2->loadSingleResidueEnergy();

        this->dt1->loadPairwiseEnergy();
        this->dt2->loadPairwiseEnergy();

    }

    void DoubleTargetDesign::mcRun(RotSequence* result1, RotSequence* result2){

         RotSequence unit1(rotGroups);
         RotSequence unit2(rotGroups.size());
         unit2.copyValue(&unit1);
         result1->copyValue(&unit1);
         result2->copyValue(&unit1);

         double energy1 = unit1.totEnergy(dt1);
         double energy2 = unit2.totEnergy(dt2);
         double energy = energy1*weight1 + energy2*weight2;

         int len = dt1->resNum;

         ResName rn;

         int randPos, randAAChoice, aaType, randChoice1, choiceIndex1, randChoice2, choiceIndex2;
         double mutEnergy;

         int aaChoiceNum, minEneChoice1, minEneChoice2;


         srand((unsigned)time(NULL));
         double minEnergy = energy;
         int count = 0;
         for(double T = T0;T>T1;T=T*annealFactor){
             bool noChange = true;
             int acNum = 0;
             for(int i=0;i<step;i++){
                 count ++;
            //     if(count %1000 == 0){
            //         printf("%-10d %10.4f\n",count,energy);
            //     }
                 randPos = rand()%len;
                 randAAChoice = rand()%posAAList.at(randPos).size();
                 aaType = posAAList.at(randPos).at(randAAChoice);

                 aaChoiceNum = posAAChoiceList[randPos][aaType].size();
                 double minMutEnergy = 99999.9;
                 for(int j=0;j<aaChoiceNum;j++){
                     mutEnergy = unit1.mutEnergy(randPos,posAAChoiceList[randPos][aaType][j], dt1);
                     if(mutEnergy < minMutEnergy){
                         minMutEnergy = mutEnergy;
                         minEneChoice1 = posAAChoiceList[randPos][aaType][j];
                     }
                 }
                 minMutEnergy = 99999.9;
                 for(int j=0;j<aaChoiceNum;j++){
                     mutEnergy = unit2.mutEnergy(randPos,posAAChoiceList[randPos][aaType][j], dt2);
                     if(mutEnergy < minMutEnergy){
                         minMutEnergy = mutEnergy;
                         minEneChoice2 = posAAChoiceList[randPos][aaType][j];
                     }
                 }

                 /*
                 choiceIndex1 = rand()%(posAAChoiceList[randPos][aaType].size());
                 choiceIndex2 = rand()%(posAAChoiceList[randPos][aaType].size());
                 randChoice1 = posAAChoiceList[randPos][aaType][choiceIndex1];
                 randChoice2 = posAAChoiceList[randPos][aaType][choiceIndex2];
                  */

                // if(count %1000 == 0)
                //     printf("pos: %3d aaChoice: %2d aaType: %c rot1: %s rot2: %s\n",randPos,randAAChoice, rn.intToSin(aaType), rotGroups.at(randPos)->rotList.at(minEneChoice1)->rotName.c_str(),  rotGroups.at(randPos)->rotList.at(minEneChoice2)->rotName.c_str());
                // cout << "pos: " << randPos << " randChoice: " << randChoice << " choiceNum: " << dt->eaList.at(randPos)->getChoiceNum() << endl;
                 mutEnergy = unit1.mutEnergy(randPos,minEneChoice1, dt1) * weight1;
                 mutEnergy += unit2.mutEnergy(randPos,minEneChoice2, dt2) * weight2;

                 if(mutEnergy < 0)
                     noChange = false;
                 if(accept(mutEnergy, T)){
                     acNum++;
                     unit1.acceptMut(randPos, minEneChoice1);
                     unit2.acceptMut(randPos, minEneChoice2);
                     energy += mutEnergy;
                     if(energy < minEnergy){
                         minEnergy = energy;
                         result1->copyValue(&unit1);
                         result2->copyValue(&unit2);
                     }
                 }
                 /*
                 if(count %10000 == 0){
                     double e1 = unit1.totEnergy(dt1);
                     double e2 = unit2.totEnergy(dt2);
                     printf("%7.4f %10.4f %10.4f %4d\n",T,e1,e2,acNum);
                 }
                 */
                // cout << "step = " << i << endl;
             }
        //     printf("%7.4f %10.4f %10.4f %4d\n",T,energy,minEnergy,acNum);
             if(noChange) break;
         }

         printf("MIN energy: %10.4f %10.4f\n",result1->totEnergy(dt1), result2->totEnergy(dt2));

    }

    void DoubleTargetDesign::doubleTargetDesignMCIdt(RotSequence* result1, RotSequence* result2, double idt, double idtFactor){

         RotSequence unit1(rotGroups);
           RotSequence unit2(rotGroups.size());


           unit2.copyValue(&unit1);

           result1->copyValue(&unit1);
           result2->copyValue(&unit1);


        double energy1 = unit1.totEnergy(dt1);
        double energy2 = unit2.totEnergy(dt2);

        int len = dt1->resNum;


        double idt12 = unit1.identity(&unit2);
        double totalEnergy = energy1 + energy2 + idtFactor*(idt12-idt)*(idt12-idt);
        double minEnergy = totalEnergy;

        printf("before mc: energy1: %10.3f energy2: %10.3f  totalEnergy: %10.3f\n",energy1, energy2, totalEnergy);

        ResName rn;

        int randPos, randAAChoice1, randAAChoice2, aaType1, aaType2;
        int randChoice1, randChoice2,  choiceIndex1, choiceIndex2;
        double mutEnergy;

        int aaChoiceNum1, aaChoiceNum2 ,minEneChoice1, minEneChoice2;

        int idtCount12;
        int oldAAType1, oldAAType2;

        double ie12Old, ie12New;

        srand((unsigned)time(NULL));

        int count = 0;
        for(double T = T0;T>T1;T=T*annealFactor){

            int acNum = 0;
            for(int i=0;i<step;i++){
                count ++;
                idtCount12 = unit1.idtCount(&unit2);
                ie12Old = idtFactor*(idtCount12*1.0/len-idt)*(idtCount12*1.0/len-idt);

                randPos = rand()%len;


                oldAAType1 =  rn.triToInt(unit1.getRotamer(randPos)->triName);
                oldAAType2 =  rn.triToInt(unit2.getRotamer(randPos)->triName);


                /*
                 * design on model1
                 */
                randAAChoice1 = rand()%posAAList[randPos].size();
                aaType1 = posAAList[randPos][randAAChoice1];

                if(oldAAType1 == oldAAType2 && aaType1 != oldAAType2){
                    ie12New = idtFactor*((idtCount12-1)*1.0/len-idt)*((idtCount12-1)*1.0/len-idt);
                }
                else if(oldAAType1 != oldAAType2 && aaType1 == oldAAType2){
                    ie12New = idtFactor*((idtCount12+1)*1.0/len-idt)*((idtCount12+1)*1.0/len-idt);
                }
                else
                    ie12New = ie12Old;


                aaChoiceNum1 = posAAChoiceList[randPos][aaType1].size();

                double minEnergy1 = 9999999.9;
                for(int k=0;k<aaChoiceNum1;k++){
                    choiceIndex1 = k;
                    if(k > 0 && rand()%100 > 50)
                        continue;
                    randChoice1 = posAAChoiceList[randPos][aaType1][choiceIndex1];
                    mutEnergy = unit1.mutEnergy(randPos, randChoice1, dt1);
                    if(mutEnergy < minEnergy1)
                    {
                        minEnergy1 = mutEnergy;
                        minEneChoice1 = randChoice1;
                    }
                }
                mutEnergy = minEnergy1 + ie12New - ie12Old;
                if(accept(mutEnergy, T)){
                    unit1.acceptMut(randPos, minEneChoice1);
                    totalEnergy += mutEnergy;
                    if(totalEnergy < minEnergy){
                        minEnergy = totalEnergy;
                        result1->copyValue(&unit1);
                    }
                }

                /*
                 * design on model2
                 */
                idtCount12 = unit1.idtCount(&unit2);
                randAAChoice2 = rand()%posAAList[randPos].size();
                aaType2 = posAAList[randPos][randAAChoice2];

                if(oldAAType2 == oldAAType1 && aaType2 != oldAAType1){
                    ie12New = idtFactor*((idtCount12-1)*1.0/len-idt)*((idtCount12-1)*1.0/len-idt);
                }
                else if(oldAAType2 != oldAAType1 && aaType2 == oldAAType1){
                    ie12New = idtFactor*((idtCount12+1)*1.0/len-idt)*((idtCount12+1)*1.0/len-idt);
                }
                else
                    ie12New = ie12Old;


                aaChoiceNum2 = posAAChoiceList[randPos][aaType2].size();
                double minEnergy2 = 9999999.9;
                for(int k=0;k<aaChoiceNum2;k++){
                    choiceIndex2 = k;
                    if(k > 0 && rand()%100 > 50)
                        continue;
                    randChoice2 = posAAChoiceList[randPos][aaType2][choiceIndex2];
                    mutEnergy = unit2.mutEnergy(randPos, randChoice2, dt2);
                    if(mutEnergy < minEnergy2)
                    {
                        minEnergy2 = mutEnergy;
                        minEneChoice2 = randChoice2;
                    }
                }
                mutEnergy = minEnergy2 + ie12New - ie12Old;
                if(accept(mutEnergy, T)){
                    unit2.acceptMut(randPos, minEneChoice2);
                    totalEnergy += mutEnergy;
                    if(totalEnergy < minEnergy){
                        minEnergy = totalEnergy;
                        result2->copyValue(&unit2);
                    }
                }




                if(count % 20000 == 0){
                    energy1 = unit1.totEnergy(dt1);
                    energy2 = unit2.totEnergy(dt2);


                    idt12 = unit1.identity(&unit2);

                    double eTot = energy1 + energy2 + idtFactor*(idt12-idt)*(idt12-idt) ;


                    printf("%7.4f %10.3f %10.3f  totEnergy: %10.3f totEnergy2: %10.3f %5.3f\n",T,energy1,energy2, totalEnergy, eTot, idt12);
                }


            }
       //     printf("%7.4f %10.4f %10.4f %4d\n",T,energy,minEnergy,acNum);

        }

        result1->copyValue(&unit1);
        result2->copyValue(&unit2);

        idt12 = result1->identity(result2);

        printf("MIN energy: %10.3f %10.3f %5.3f\n",result1->totEnergy(dt1), result2->totEnergy(dt2), idt12);


    }

    void DoubleTargetDesign::threeTargetDesignMC(RotSequence* result1, RotSequence* result2, RotSequence* result3, RotSequence* result4, double idt, double idtFactor){
            /*
             * unit1: sequence design based on model 1
             * unit2: sequence design based on model 2
             * unit3: sequence design based on both models, rotamers are built on model1
             * unit4: sequence design based on both models, rotamers are built on model2
             * amino acid sequences of unit3 and unit4 are the same
             * amino acid sequences identity between unit1 and unit3 is around idt
             * amino acid sequences identity between unit2 and unit3 is around idt
             */

         RotSequence unit1(rotGroups);
            RotSequence unit2(rotGroups.size());
         RotSequence unit3(rotGroups.size());
         RotSequence unit4(rotGroups.size());

            unit2.copyValue(&unit1);
            unit3.copyValue(&unit1);
            unit4.copyValue(&unit1);
            result1->copyValue(&unit1);
            result2->copyValue(&unit1);
            result3->copyValue(&unit1);
            result4->copyValue(&unit1);

         double energy1 = unit1.totEnergy(dt1);
         double energy2 = unit2.totEnergy(dt2);
         double energy3 = energy1;
         double energy4 = energy2;

         int len = dt1->resNum;


         double idt13 = unit1.identity(&unit3);
         double idt23 = unit2.identity(&unit3);
         double totalEnergy = energy1 + energy2 + 0.5*(energy3+energy4) + idtFactor*(idt13-idt)*(idt13-idt) + idtFactor*(idt23-idt)*(idt23-idt);
         double minEnergy = totalEnergy;

         printf("before mc: energy1: %10.3f energy2: %10.3f energy3: %10.3f energy4: %10.3f totalEnergy: %10.3f\n",energy1, energy2, energy3, energy4, totalEnergy);

         ResName rn;

         int randPos, randAAChoice1, randAAChoice2, randAAChoice34, aaType1, aaType2, aaType34;
         int randChoice1, randChoice2, randChoice3, randChoice4, choiceIndex1, choiceIndex2, choiceIndex3, choiceIndex4;
         double mutEnergy;

         int aaChoiceNum1, aaChoiceNum2, aaChoiceNum34,minEneChoice1, minEneChoice2, minEneChoice3, minEneChoice4;

         int idtCount13, idtCount23;
         int oldAAType1, oldAAType2, oldAAType34;

         double ie13Old, ie23Old, ie13New, ie23New;

         srand((unsigned)time(NULL));

         int count = 0;
         for(double T = T0;T>T1;T=T*annealFactor){

             int acNum = 0;
             for(int i=0;i<step;i++){
                 count ++;
                 idtCount13 = unit1.idtCount(&unit3);
                 idtCount23 = unit2.idtCount(&unit3);
                 ie13Old = idtFactor*(idtCount13*1.0/len-idt)*(idtCount13*1.0/len-idt);
                 ie23Old = idtFactor*(idtCount23*1.0/len-idt)*(idtCount23*1.0/len-idt);

                 randPos = rand()%len;


                 oldAAType1 =  rn.triToInt(unit1.getRotamer(randPos)->triName);
                 oldAAType2 =  rn.triToInt(unit2.getRotamer(randPos)->triName);
                 oldAAType34 =  rn.triToInt(unit3.getRotamer(randPos)->triName);


                 /*
                  * design on model1
                  */
                 randAAChoice1 = rand()%posAAList[randPos].size();
                 aaType1 = posAAList[randPos][randAAChoice1];

                 if(oldAAType1 == oldAAType34 && aaType1 != oldAAType1){
                     ie13New = idtFactor*((idtCount13-1)*1.0/len-idt)*((idtCount13-1)*1.0/len-idt);
                 }
                 else if(oldAAType1 != oldAAType34 && aaType1 == oldAAType34){
                     ie13New = idtFactor*((idtCount13+1)*1.0/len-idt)*((idtCount13+1)*1.0/len-idt);
                 }
                 else
                     ie13New = ie13Old;


                 aaChoiceNum1 = posAAChoiceList[randPos][aaType1].size();

                 double minEnergy1 = 9999999.9;
                 for(int k=0;k<aaChoiceNum1;k++){
                     choiceIndex1 = k;
                     if(k > 0 && rand()%100 > 50)
                         continue;
                     randChoice1 = posAAChoiceList[randPos][aaType1][choiceIndex1];
                     mutEnergy = unit1.mutEnergy(randPos, randChoice1, dt1);
                     if(mutEnergy < minEnergy1)
                     {
                         minEnergy1 = mutEnergy;
                         minEneChoice1 = randChoice1;
                     }
                 }
                 mutEnergy = minEnergy1 + ie13New - ie13Old;
                 if(accept(mutEnergy, T)){
                     unit1.acceptMut(randPos, minEneChoice1);
                     totalEnergy += mutEnergy;
                     if(totalEnergy < minEnergy){
                         minEnergy = totalEnergy;
                         result1->copyValue(&unit1);
                     }
                 }

                 /*
                  * design on model2
                  */
                 randAAChoice2 = rand()%posAAList[randPos].size();
                 aaType2 = posAAList[randPos][randAAChoice2];

                 if(oldAAType2 == oldAAType34 && aaType2 != oldAAType2){
                     ie23New = idtFactor*((idtCount23-1)*1.0/len-idt)*((idtCount23-1)*1.0/len-idt);
                 }
                 else if(oldAAType2 != oldAAType34 && aaType2 == oldAAType34){
                     ie23New = idtFactor*((idtCount23+1)*1.0/len-idt)*((idtCount23+1)*1.0/len-idt);
                 }
                 else
                     ie23New = ie23Old;


                 aaChoiceNum2 = posAAChoiceList[randPos][aaType2].size();
                 double minEnergy2 = 9999999.9;
                 for(int k=0;k<aaChoiceNum2;k++){
                     choiceIndex2 = k;
                     if(k > 0 && rand()%100 > 50)
                         continue;
                     randChoice2 = posAAChoiceList[randPos][aaType2][choiceIndex2];
                     mutEnergy = unit2.mutEnergy(randPos, randChoice2, dt2);
                     if(mutEnergy < minEnergy2)
                     {
                         minEnergy2 = mutEnergy;
                         minEneChoice2 = randChoice2;
                     }
                 }
                 mutEnergy = minEnergy2 + ie23New - ie23Old;
                 if(accept(mutEnergy, T)){
                     unit2.acceptMut(randPos, minEneChoice2);
                     totalEnergy += mutEnergy;
                     if(totalEnergy < minEnergy){
                         minEnergy = totalEnergy;
                         result2->copyValue(&unit2);
                     }
                 }


                 /*
                  * design on model34
                  */
                 idtCount13 = unit1.idtCount(&unit3);
                 idtCount23 = unit2.idtCount(&unit3);
                 ie13Old = idtFactor*(idtCount13*1.0/len-idt)*(idtCount13*1.0/len-idt);
                 ie23Old = idtFactor*(idtCount23*1.0/len-idt)*(idtCount23*1.0/len-idt);

                 oldAAType1 =  rn.triToInt(unit1.getRotamer(randPos)->triName);
                 oldAAType2 =  rn.triToInt(unit2.getRotamer(randPos)->triName);
                 oldAAType34 =  rn.triToInt(unit3.getRotamer(randPos)->triName);

                 randAAChoice34 = rand()%posAAList[randPos].size();
                 aaType34 = posAAList[randPos][randAAChoice34];

                 ie13New = ie13Old;
                 ie23New = ie23Old;
                 if(oldAAType34 != oldAAType1 && aaType34 == oldAAType1)
                     ie13New = idtFactor*((idtCount13+1)*1.0/len-idt)*((idtCount13+1)*1.0/len-idt);
                 if(oldAAType34 == oldAAType1 && aaType34 != oldAAType1)
                     ie13New = idtFactor*((idtCount13-1)*1.0/len-idt)*((idtCount13-1)*1.0/len-idt);
                 if(oldAAType34 != oldAAType2 && aaType34 == oldAAType2)
                     ie23New = idtFactor*((idtCount23+1)*1.0/len-idt)*((idtCount23+1)*1.0/len-idt);
                 if(oldAAType34 == oldAAType2 && aaType34 != oldAAType2)
                     ie23New = idtFactor*((idtCount23-1)*1.0/len-idt)*((idtCount23-1)*1.0/len-idt);

                 aaChoiceNum34 = posAAChoiceList[randPos][aaType34].size();
                 double minMutEnergy = 99999.9;
                 for(int j=0;j<aaChoiceNum34;j++){
                     if(j > 0 && rand()%100 > 70)
                         continue;
                     mutEnergy = unit3.mutEnergy(randPos,posAAChoiceList[randPos][aaType34][j], dt1);
                     if(mutEnergy < minMutEnergy){
                         minMutEnergy = mutEnergy;
                         minEneChoice3 = posAAChoiceList[randPos][aaType34][j];
                     }
                 }
                 minMutEnergy = 99999.9;
                 for(int j=0;j<aaChoiceNum34;j++){
                     if(j > 0 && rand()%100 > 70)
                         continue;
                     mutEnergy = unit4.mutEnergy(randPos,posAAChoiceList[randPos][aaType34][j], dt2);
                     if(mutEnergy < minMutEnergy){
                         minMutEnergy = mutEnergy;
                         minEneChoice4 = posAAChoiceList[randPos][aaType34][j];
                     }
                 }

                 /*
                 choiceIndex1 = rand()%(posAAChoiceList[randPos][aaType].size());
                 choiceIndex2 = rand()%(posAAChoiceList[randPos][aaType].size());
                 randChoice1 = posAAChoiceList[randPos][aaType][choiceIndex1];
                 randChoice2 = posAAChoiceList[randPos][aaType][choiceIndex2];
                  */

                // if(count %1000 == 0)
                //     printf("pos: %3d aaChoice: %2d aaType: %c rot1: %s rot2: %s\n",randPos,randAAChoice, rn.intToSin(aaType), rotGroups.at(randPos)->rotList.at(minEneChoice1)->rotName.c_str(),  rotGroups.at(randPos)->rotList.at(minEneChoice2)->rotName.c_str());
                // cout << "pos: " << randPos << " randChoice: " << randChoice << " choiceNum: " << dt->eaList.at(randPos)->getChoiceNum() << endl;
                 mutEnergy = 0.5*unit3.mutEnergy(randPos,minEneChoice3, dt1) * weight1;
                 mutEnergy += 0.5*unit4.mutEnergy(randPos,minEneChoice4, dt2) * weight2;
                 mutEnergy += ie13New - ie13Old + ie23New - ie23Old;


                 if(accept(mutEnergy, T)){
                     acNum++;
                     unit3.acceptMut(randPos, minEneChoice3);
                     unit4.acceptMut(randPos, minEneChoice4);
                     totalEnergy += mutEnergy;
                     if(totalEnergy < minEnergy){
                         minEnergy = totalEnergy;
                         result3->copyValue(&unit3);
                         result4->copyValue(&unit4);
                     }
                 }
                 /*
                 if(count % 10000 == 0){
                     energy1 = unit1.totEnergy(dt1);
                     energy2 = unit2.totEnergy(dt2);
                     energy3 = unit3.totEnergy(dt1);
                     energy4 = unit4.totEnergy(dt2);

                     idt13 = unit1.identity(&unit3);
                     idt23 = unit2.identity(&unit3);
                     double eTot = energy1 + energy2 + 0.5*(energy3+energy4) + idtFactor*(idt13-idt)*(idt13-idt) + idtFactor*(idt23-idt)*(idt23-idt);


                     printf("%7.4f %10.3f %10.3f %10.3f %10.3f totEnergy: %10.3f totEnergy2: %10.3f %5.3f %5.3f\n",T,energy1,energy2,energy3, energy4, totalEnergy, eTot, idt13, idt23);
                 }
                 */
                // cout << "step = " << i << endl;
             }
        //     printf("%7.4f %10.4f %10.4f %4d\n",T,energy,minEnergy,acNum);

         }

         result1->copyValue(&unit1);
         result2->copyValue(&unit2);
         result3->copyValue(&unit3);
         result4->copyValue(&unit4);

         idt13 = result1->identity(result3);
         idt23 = result2->identity(result3);

         printf("MIN energy: %10.3f %10.3f %10.3f %10.3f %5.3f %5.3f\n",result1->totEnergy(dt1), result2->totEnergy(dt2), result3->totEnergy(dt1), result4->totEnergy(dt2), idt13, idt23);



    }

    void DoubleTargetDesign::threeTargetDesignMC34Fixed(RotSequence* result1, RotSequence* result2, RotSequence* result3, RotSequence* result4, double idt, double idtFactor){
            /*
             * unit1: sequence design based on model 1
             * unit2: sequence design based on model 2
             * unit3: sequence design based on both models, rotamers are built on model1
             * unit4: sequence design based on both models, rotamers are built on model2
             * amino acid sequences of unit3 and unit4 are the same
             * amino acid sequences identity between unit1 and unit3 is around idt
             * amino acid sequences identity between unit2 and unit3 is around idt
             */

         RotSequence unit1(rotGroups.size());
            RotSequence unit2(rotGroups.size());
         RotSequence unit3(rotGroups.size());
         RotSequence unit4(rotGroups.size());

         unit1.copyValue(result3);
            unit2.copyValue(&unit1);
            unit3.copyValue(result3);
            unit4.copyValue(result4);

            result1->copyValue(&unit1);
            result2->copyValue(&unit1);


         double energy1 = unit1.totEnergy(dt1);
         double energy2 = unit2.totEnergy(dt2);
         double energy3 = energy1;
         double energy4 = energy2;

         int len = dt1->resNum;


         double idt13 = unit1.identity(&unit3);
         double idt23 = unit2.identity(&unit3);
         double totalEnergy = energy1 + energy2 + 0.5*(energy3+energy4) + idtFactor*(idt13-idt)*(idt13-idt) + idtFactor*(idt23-idt)*(idt23-idt);
         double minEnergy = totalEnergy;

         printf("before mc: energy1: %10.3f energy2: %10.3f energy3: %10.3f energy4: %10.3f totalEnergy: %10.3f\n",energy1, energy2, energy3, energy4, totalEnergy);

         ResName rn;

         int randPos, randAAChoice1, randAAChoice2, randAAChoice34, aaType1, aaType2, aaType34;
         int randChoice1, randChoice2, randChoice3, randChoice4, choiceIndex1, choiceIndex2, choiceIndex3, choiceIndex4;
         double mutEnergy;

         int aaChoiceNum1, aaChoiceNum2, aaChoiceNum34,minEneChoice1, minEneChoice2, minEneChoice3, minEneChoice4;

         int idtCount13, idtCount23;
         int oldAAType1, oldAAType2, oldAAType34;

         double ie13Old, ie23Old, ie13New, ie23New;

         srand((unsigned)time(NULL));

         int count = 0;
         for(double T = 1.0;T>0.001;T=T*annealFactor){
             int acNum = 0;
             for(int i=0;i<step;i++){
                 count ++;
                 idtCount13 = unit1.idtCount(&unit3);
                 idtCount23 = unit2.idtCount(&unit3);
                 ie13Old = idtFactor*(idtCount13*1.0/len-idt)*(idtCount13*1.0/len-idt);
                 ie23Old = idtFactor*(idtCount23*1.0/len-idt)*(idtCount23*1.0/len-idt);
                 randPos = rand()%len;
                 oldAAType1 =  rn.triToInt(unit1.getRotamer(randPos)->triName);
                 oldAAType2 =  rn.triToInt(unit2.getRotamer(randPos)->triName);
                 oldAAType34 =  rn.triToInt(unit3.getRotamer(randPos)->triName);


                 /*
                  * design on model1
                  */
                 randAAChoice1 = rand()%posAAList[randPos].size();
                 aaType1 = posAAList[randPos][randAAChoice1];

                 if(oldAAType1 == oldAAType34 && aaType1 != oldAAType1){
                     ie13New = idtFactor*((idtCount13-1)*1.0/len-idt)*((idtCount13-1)*1.0/len-idt);
                 }
                 else if(oldAAType1 != oldAAType34 && aaType1 == oldAAType34){
                     ie13New = idtFactor*((idtCount13+1)*1.0/len-idt)*((idtCount13+1)*1.0/len-idt);
                 }
                 else
                     ie13New = ie13Old;


                 aaChoiceNum1 = posAAChoiceList[randPos][aaType1].size();

                 double minEnergy1 = 9999999.9;
                 for(int k=0;k<aaChoiceNum1;k++){
                     choiceIndex1 = k;
                     if(k > 0 && rand()%100 > 50)
                         continue;
                     randChoice1 = posAAChoiceList[randPos][aaType1][choiceIndex1];
                     mutEnergy = unit1.mutEnergy(randPos, randChoice1, dt1);
                     if(mutEnergy < minEnergy1)
                     {
                         minEnergy1 = mutEnergy;
                         minEneChoice1 = randChoice1;
                     }
                 }
                 mutEnergy = minEnergy1 + ie13New - ie13Old;
                 if(accept(mutEnergy, T)){
                     unit1.acceptMut(randPos, minEneChoice1);
                     totalEnergy += mutEnergy;
                     if(totalEnergy < minEnergy){
                         minEnergy = totalEnergy;
                         result1->copyValue(&unit1);
                     }
                 }

                 /*
                  * design on model2
                  */
                 randAAChoice2 = rand()%posAAList[randPos].size();
                 aaType2 = posAAList[randPos][randAAChoice2];

                 if(oldAAType2 == oldAAType34 && aaType2 != oldAAType2){
                     ie23New = idtFactor*((idtCount23-1)*1.0/len-idt)*((idtCount23-1)*1.0/len-idt);
                 }
                 else if(oldAAType2 != oldAAType34 && aaType2 == oldAAType34){
                     ie23New = idtFactor*((idtCount23+1)*1.0/len-idt)*((idtCount23+1)*1.0/len-idt);
                 }
                 else
                     ie23New = ie23Old;

                 aaChoiceNum2 = posAAChoiceList[randPos][aaType2].size();
                 double minEnergy2 = 9999999.9;
                 for(int k=0;k<aaChoiceNum2;k++){
                     choiceIndex2 = k;
                     if(k > 0 && rand()%100 > 50)
                         continue;
                     randChoice2 = posAAChoiceList[randPos][aaType2][choiceIndex2];
                     mutEnergy = unit2.mutEnergy(randPos, randChoice2, dt2);
                     if(mutEnergy < minEnergy2)
                     {
                         minEnergy2 = mutEnergy;
                         minEneChoice2 = randChoice2;
                     }
                 }
                 mutEnergy = minEnergy2 + ie23New - ie23Old;
                 if(accept(mutEnergy, T)){
                     unit2.acceptMut(randPos, minEneChoice2);
                     totalEnergy += mutEnergy;
                     if(totalEnergy < minEnergy){
                         minEnergy = totalEnergy;
                         result2->copyValue(&unit2);
                     }
                 }

                 /*
                 if(count % 10000 == 0){
                     energy1 = unit1.totEnergy(dt1);
                     energy2 = unit2.totEnergy(dt2);
                     energy3 = unit3.totEnergy(dt1);
                     energy4 = unit4.totEnergy(dt2);

                     idt13 = unit1.identity(&unit3);
                     idt23 = unit2.identity(&unit3);
                     double eTot = energy1 + energy2 + 0.5*(energy3+energy4) + idtFactor*(idt13-idt)*(idt13-idt) + idtFactor*(idt23-idt)*(idt23-idt);


                     printf("%7.4f %10.3f %10.3f %10.3f %10.3f totEnergy: %10.3f totEnergy2: %10.3f %5.3f %5.3f\n",T,energy1,energy2,energy3, energy4, totalEnergy, eTot, idt13, idt23);
                 }
                 */
             }

         }

         result1->copyValue(&unit1);
         result2->copyValue(&unit2);
         result3->copyValue(&unit3);
         result4->copyValue(&unit4);


         idt13 = result1->identity(result3);
         idt23 = result2->identity(result3);

         printf("MIN energy: %10.3f %10.3f %10.3f %10.3f %5.3f %5.3f\n",result1->totEnergy(dt1), result2->totEnergy(dt2), result3->totEnergy(dt1), result4->totEnergy(dt2), idt13, idt23);


    }

    void DoubleTargetDesign::printPDB(RotSequence* unit, string outputFile, DesignTemplate* dt){
        ProteinChain pc;
        char s[20];
        for(int i=0;i<dt->bsList.size();i++){
            BackBoneSite* bs = dt->bsList.at(i);
            Rotamer* rot = unit->getRotamer(i);
            Residue* res = new Residue(bs->strresid, bs->chainid, rot->triName);
            res->addAtom(new Atom("N", bs->ncrd()));
            res->addAtom(new Atom("CA", bs->cacrd()));
            res->addAtom(new Atom("C", bs->ccrd()));
            res->addAtom(new Atom("O", bs->ocrd()));
            res->buildRotamer(rot);
            pc.addResidue(res);
        }
        ofstream output(outputFile, ios::out);
        if(!output.is_open())
        {
            cout << "fail to open file " << outputFile << endl;
            exit(1);
        }
        pc.printPDBFormat(output, 1);

        Residue* p1;
        Atom* p2;
        for(int i=0;i<pc.getChainLength();i++){
            p1 = pc.getResList().at(i);
            for(int j=0;j<p1->getAtomList()->size();j++){
                p2 = p1->getAtomList()->at(j);
                delete p2;
            }
            delete p1;
        }
        output.close();
    }

    DoubleTargetDesign::~DoubleTargetDesign(){
        for(int i=0;i<rotGroups.size();i++){
            delete rotGroups[i];
        }
    }



    /*
     * added by xuyang, 2019.3.31
     */

    RotSequence::RotSequence(vector<RotamerGroup*>& groups, std::vector<AANumberControlInMC> &aacmc){

        std::cout <<"Starting Initiating Sequence..." <<std::endl;

        int len = groups.size();
        this->seqLength = len;
        this->seqChoices = new int[len];
        this->aaSeq = new char[len];
        for(int i=0;i<len;i++){
            this->rotGroups.push_back(groups[i]);
            rotGroups.back()->updataidinlist();
        }

        std::cout <<"Check User-Defined Amino Acid Proprotion..." <<std::endl;

        for(const auto &ac:aacmc) {
            int nt=0;
            for(int i:ac.sites) {
                if(rotGroups[i]->containaa(ac.aa)) nt++;
            }
            if(nt>=ac.min) continue;
            std::cout <<"Can Not Find Solution Satisfied Count of Amino Acid!\t\t" <<ac.aa <<std::endl;
            exit(1);
        }

        std::cout <<"Make Initial Random Sequence..." <<std::endl;

        setRandomChoice();

        bool isgood{false};
        for(int i=0;i<10000;i++) {
            std::cout <<"Adjust Sequence " <<i <<std::endl;
            bool allissatisfied{true};
            for(int j=0;j<aacmc.size();j++) {
                if(adjustSeq_1time(aacmc[j])) {
                    allissatisfied = false;
                    break;
                }
            }
            if(allissatisfied) {
                isgood=true;
                break;
            }
        }

        if(!isgood) {
            std::cout <<"Can Not Find Solution Satisfied Count of Amino Acid!" <<std::endl;
            exit(1);
        } else {
            std::cout <<"Initial Sequence Is:" <<std::endl;
            std::string seq=aaSeq;
            std::cout <<'\t' <<seq <<std::endl;
        }
    }

    bool RotSequence::adjustSeq_1time(AANumberControlInMC &ac) {
        int nt=0;
        for(int i:ac.sites) {
            if(aaSeq[i]==ac.aa) nt++;
        }
        ac.num = nt;
        std::cout <<ac.num <<'\t' <<ac.min <<'\t' <<ac.max <<std::endl;
        //if(nt<ac.min) return '-';
        //else if(nt>ac.max) return '+';
        //else return '0';
        if(nt<ac.min) {
            std::vector<int> vi;
            for(int i:ac.sites) if(aaSeq[i]!=ac.aa) if(rotGroups[i]->containaa(ac.aa)) vi.push_back(i);
            if(vi.empty()) {
                std::cout <<"Min AA Proportion will change:  " <<ac.aa <<"    " <<ac.min <<"  ->  " <<nt <<std::endl;
                for(auto &s:ac.sites) std::cout <<'\t' <<s;
                std::cout <<endl;
                ac.min = nt;
                return true;
            }
            int i1=rand()%(vi.size());
            int pos=vi[i1];
            std::vector<int> vi1=goodChoiceOnly(ac.aa,pos);
            assert(vi1.size()>0);
            int i2=rand()%(vi1.size());
            int randchose = vi1[i2];

            this->seqChoices[pos] = randchose;
            Rotamer* rot = this->rotGroups[pos]->rotList[this->seqChoices[pos]];
            aaSeq[pos] = rn.triToSin(rot->triName);

            ac.num ++;

            return true;
        } else if(nt>ac.max) {
            std::vector<int> vi;
            for(int i:ac.sites) if(aaSeq[i]==ac.aa) if(rotGroups[i]->containother(ac.aa)) vi.push_back(i);
            if(vi.empty()) {
                //std::cout <<"Can Not Adjust Sequence!" <<std::endl;
                //exit(1);
                std::cout <<"Max AA Proportion will change:  " <<ac.aa <<"    " <<ac.max <<"  ->  " <<nt <<std::endl;
                for(auto &s:ac.sites) std::cout <<'\t' <<s;
                std::cout <<endl;
                ac.max = nt;
                return true;
            }
            int i1=rand()%(vi.size());
            int pos=vi[i1];
            std::vector<int> vi1=goodChoiceExcept(ac.aa,pos);
            assert(vi1.size()>0);
            int i2=rand()%(vi1.size());
            int randchose = vi1[i2];

            this->seqChoices[pos] = randchose;
            Rotamer* rot = this->rotGroups[pos]->rotList[this->seqChoices[pos]];
            aaSeq[pos] = rn.triToSin(rot->triName);

            ac.num --;

            return true;
        } else {
            return false;
        }
    }

    std::vector<int> RotSequence::goodChoiceOnly(char aa, int pos) {
        std::vector<int> vi;
        for(int i=0;i<rotGroups[pos]->rotList.size();i++) {
            if(rn.triToSin(rotGroups[pos]->rotList[i]->triName)==aa) vi.push_back(i);
        }
        return vi;
    }
    std::vector<int> RotSequence::goodChoiceExcept(char aa, int pos) {
        std::vector<int> vi;
        for(int i=0;i<rotGroups[pos]->rotList.size();i++) {
            if(rn.triToSin(rotGroups[pos]->rotList[i]->triName)!=aa) vi.push_back(i);
        }
        return vi;
    }

    void RotSequence::calcEnergyComponents(const DesignTemplate* dt, double* s1,
                                           double* rot, double* ref) const {
        assert(dt->eaList.size() == seqLength);
        *s1 = 0;
        *rot = 0;
        *ref = 0;
        for (int i = 0; i < seqLength; i++) {
            *s1 += dt->ecList.at(i).s1.at(seqChoices[i]);
            *rot += dt->ecList.at(i).rot.at(seqChoices[i]);
            *ref += dt->ecList.at(i).ref.at(seqChoices[i]);
        }
    }

    double DesignMC::mcRunWithAACountRestrain(RotSequence* result,
            std::vector<AANumberControlInMC> &aacmc){

         RotSequence unit(dt->rotGroups, aacmc);
         result->copyValue(&unit);
         double energy = unit.totEnergy(dt);
         int len = unit.getLength();
         int randPos, randChoice;
         double mutEnergy;

         std::map<int,std::vector<AANumberControlInMC*>> aanc;
         for(int i=0;i<aacmc.size();i++) {
             for(int j:aacmc[i].sites) {
                 if(aanc.find(j)==aanc.end()) aanc.insert({j,std::vector<AANumberControlInMC*>()});
                 aanc.at(j).push_back(&(aacmc[i]));
             }
         }

         srand((unsigned)time(NULL));
         double minEnergy = energy;
         int count = 0;
         for(double T = T0;T>T1;T=T*annealFactor){
             bool noChange = true;
             int acNum = 0;
             //std::cout <<"Temperature: " <<T <<std::endl;
             for(int i=0;i<step;i++){
                 count ++;
                 randPos = rand()%len;

                 //randChoice = rand()%(unit.choiceNum(randPos));

                 char aa_old;
                 char aa_new;

                 if(aanc.find(randPos)==aanc.end()) randChoice = rand()%(unit.choiceNum(randPos));
                 else {
                     std::set<int> tobedelete;
                     aa_old=unit.getaatype(randPos);
                     RotamerGroup* rgs=unit.getRotGroup(randPos);
                     for(AANumberControlInMC* &ac:aanc.at(randPos)) {
                         if(ac->min==ac->max) {
                               for(int j=0;j<rgs->rotList.size();j++) {
                                     char c=unit.getrn().triToSin(rgs->rotList[j]->triName);
                                     if(c!=aa_old) tobedelete.insert(j);
                                }
                         } else if(aa_old==ac->aa) {
                             if(ac->min==ac->num) {
                                    for(int j=0;j<rgs->rotList.size();j++) {
                                         char c=unit.getrn().triToSin(rgs->rotList[j]->triName);
                                         if(c!=aa_old) tobedelete.insert(j);
                                    }
                             }
                         } else {
                             if(ac->max==ac->num) {
                                 for(int j=0;j<rgs->rotList.size();j++) {
                                     char c=unit.getrn().triToSin(rgs->rotList[j]->triName);
                                     if(c==ac->aa) tobedelete.insert(j);
                                 }
                             }
                         }
                     }

                     int nt=unit.choiceNum(randPos);
                     std::vector<int> tobechose;
                     for(int j=0;j<nt;j++) if(tobedelete.find(j)==tobedelete.end()) tobechose.push_back(j);
                     assert(tobechose.size()>0);
                     int ix = rand()%(tobechose.size());
                     randChoice = tobechose[ix];

                     aa_new = unit.getrn().triToSin(rgs->rotList[randChoice]->triName);
                 }

                 mutEnergy = unit.mutEnergy(randPos,randChoice, dt);
                 if(mutEnergy < 0)
                     noChange = false;
                 if(accept(mutEnergy, T)){
                     acNum++;
                     unit.acceptMutAndUpdateSeq(randPos, randChoice);
                     energy += mutEnergy;
                     if(energy < minEnergy){
                         minEnergy = energy;
                         result->copyValue(&unit);
                     }
                     if(aa_old!=aa_new && aanc.find(randPos)!=aanc.end()) {
                         for(int j=0;j<aanc.at(randPos).size();j++) {
                             if(aanc.at(randPos)[j]->aa==aa_old) aanc.at(randPos)[j]->num--;
                             if(aanc.at(randPos)[j]->aa==aa_new) aanc.at(randPos)[j]->num++;
                         }
                     }
                 }
             }
             printf("%7.4f %10.4f %10.4f %4d\n",T,energy,minEnergy,acNum);
             if(noChange) break;
         }
         double finalEnergy = result->totEnergy(dt);

         return finalEnergy;
    }

    double DesignMC::mcRunWithIdtAACountRestrain(RotSequence* result, vector<string>& preSeqs,
            double idtCutoff,double sigma, std::vector<AANumberControlInMC> &aacmc){

         RotSequence unit(dt->rotGroups, aacmc);
         result->copyValue(&unit);
         double energy = unit.totEnergy(dt);
         int len = unit.getLength();
         int randPos, randChoice;
         double mutEnergy;
         double idtMutEnergy;

         std::map<int,std::vector<AANumberControlInMC*>> aanc;
         for(int i=0;i<aacmc.size();i++) {
             for(int j:aacmc[i].sites) {
                 if(aanc.find(j)==aanc.end()) aanc.insert({j,std::vector<AANumberControlInMC*>()});
                 aanc.at(j).push_back(&(aacmc[i]));
             }
         }

         srand((unsigned)time(NULL));
         double minEnergy = energy;
         int count = 0;
         for(double T = T0;T>T1;T=T*annealFactor){
             bool noChange = true;
             int acNum = 0;
             //std::cout <<"Temperature: " <<T <<std::endl;
             for(int i=0;i<step;i++){
                 count ++;
                 randPos = rand()%len;

                 //randChoice = rand()%(unit.choiceNum(randPos));

                 char aa_old;
                 char aa_new;

                 if(aanc.find(randPos)==aanc.end()) randChoice = rand()%(unit.choiceNum(randPos));
                 else {
                     std::set<int> tobedelete;
                     aa_old=unit.getaatype(randPos);
                     RotamerGroup* rgs=unit.getRotGroup(randPos);
                     for(AANumberControlInMC* &ac:aanc.at(randPos)) {
                         if(ac->min==ac->max) {
                               for(int j=0;j<rgs->rotList.size();j++) {
                                     char c=unit.getrn().triToSin(rgs->rotList[j]->triName);
                                     if(c!=aa_old) tobedelete.insert(j);
                                }
                         } else if(aa_old==ac->aa) {
                             if(ac->min==ac->num) {
                                    for(int j=0;j<rgs->rotList.size();j++) {
                                         char c=unit.getrn().triToSin(rgs->rotList[j]->triName);
                                         if(c!=aa_old) tobedelete.insert(j);
                                    }
                             }
                         } else {
                             if(ac->max==ac->num) {
                                 for(int j=0;j<rgs->rotList.size();j++) {
                                     char c=unit.getrn().triToSin(rgs->rotList[j]->triName);
                                     if(c==ac->aa) tobedelete.insert(j);
                                 }
                             }
                         }
                     }

                     int nt=unit.choiceNum(randPos);
                     std::vector<int> tobechose;
                     for(int j=0;j<nt;j++) if(tobedelete.find(j)==tobedelete.end()) tobechose.push_back(j);
                     assert(tobechose.size()>0);
                     int ix = rand()%(tobechose.size());
                     randChoice = tobechose[ix];

                     aa_new = unit.getrn().triToSin(rgs->rotList[randChoice]->triName);
                 }

                 mutEnergy = unit.mutEnergy(randPos,randChoice, dt);
                 idtMutEnergy = unit.mutIdentityEnergy(preSeqs,randPos,randChoice,idtCutoff,sigma);
                 mutEnergy += idtMutEnergy;
                 if(mutEnergy < 0)
                     noChange = false;
                 if(accept(mutEnergy, T)){
                     acNum++;
                     unit.acceptMutAndUpdateSeq(randPos, randChoice);
                     energy += mutEnergy;
                     if(energy < minEnergy){
                         minEnergy = energy;
                         result->copyValue(&unit);
                     }
                     if(aa_old!=aa_new && aanc.find(randPos)!=aanc.end()) {
                         for(int j=0;j<aanc.at(randPos).size();j++) {
                             if(aanc.at(randPos)[j]->aa==aa_old) aanc.at(randPos)[j]->num--;
                             if(aanc.at(randPos)[j]->aa==aa_new) aanc.at(randPos)[j]->num++;
                         }
                     }
                 }
             }
             printf("%7.4f %10.4f %10.4f %4d\n",T,energy,minEnergy,acNum);
             if(noChange) break;
         }
         double finalEnergy = result->totEnergy(dt);

         return finalEnergy;
    }


    /*
    bool RotSequence::AACountRestrain(const std::vector<AANumberControlInMC> &aacmc,
            const std::map<char,std::vector<std::pair<std::vector<int>,std::set<int>>>> &aacr){
        std::set<char> a20{'G','A','V','L','I','S','T','C','M','D','E','N','Q','K','R','H','P','F','Y','W'};
        assert(a20.size()==20);
        for(const AANumberControlInMC &ac:aacmc) {

        }
        for(const auto &ac:aacr) {
            for(int i=0;i<ac.second.size();i++) {
                const std::set<int>& as = ac.second[i].second;
                int nt=0;
                for(auto & i: as) {
                    if(aaSeq[i]==ac.first) nt++;
                }
                int min=ac.second[i].first[0];
                int max=ac.second[i].first[1];
                if(max<min) {
                    min=ac.second[i].first[1];
                    max=ac.second[i].first[0];
                }
                //std::cout <<nt <<'\t' <<min <<'\t' <<max <<std::endl;
                if(nt<min) return false;
                if(nt>max) return false;
            }
        }
        return true;
    }

    bool RotSequence::AACountRestrain(
            const std::map<char,std::vector<std::pair<std::vector<int>,std::set<int>>>> &aacr){
        for(const auto &ac:aacr) {
            for(int i=0;i<ac.second.size();i++) {
                const std::set<int>& as = ac.second[i].second;
                int nt=0;
                for(auto & i: as) {
                    if(aaSeq[i]==ac.first) nt++;
                }
                int min=ac.second[i].first[0];
                int max=ac.second[i].first[1];
                if(max<min) {
                    min=ac.second[i].first[1];
                    max=ac.second[i].first[0];
                }
                //std::cout <<nt <<'\t' <<min <<'\t' <<max <<std::endl;
                if(nt<min) return false;
                if(nt>max) return false;
            }
        }
        return true;
    }

    RotSequence::RotSequence(vector<RotamerGroup*>& groups,
            const std::map<char,std::vector<std::pair<std::vector<int>,std::set<int>>>> &aacr){
        int len = groups.size();
        this->seqLength = len;
        this->seqChoices = new int[len];
        this->aaSeq = new char[len];
        for(int i=0;i<len;i++){
            this->rotGroups.push_back(groups[i]);
        }
        bool isgood{false};
        for(int i=0;i<1000;i++) {
            setRandomChoice();
            if(!AACountRestrain(aacr)) continue;
            isgood=true;
            break;
        }
        if(!isgood) {
            std::cout <<"Can Not Find Solution Satisfied Count of Amino Acid!" <<std::endl;
            exit(1);
        }
    }

    bool RotSequence::AACountRestrain(int randPos, int randChoice,
            const std::map<char,std::vector<std::pair<std::vector<int>,std::set<int>>>> &aacr){
        char oldaa = aaSeq[randPos];
        Rotamer* rot = this->rotGroups[randPos]->rotList[randChoice];
        char newaa = rn.triToSin(rot->triName);
        if(oldaa==newaa) return true;
        for(const auto &ac:aacr) {
            if(ac.first!=oldaa && ac.first!=newaa) continue;
            for(int i=0;i<ac.second.size();i++) {
                const std::set<int>& as = ac.second[i].second;
                int nt=0;
                for(auto & i: as) {
                    if(aaSeq[i]==ac.first) nt++;
                }
                int min=ac.second[i].first[0];
                int max=ac.second[i].first[1];
                if(max<min) {
                    min=ac.second[i].first[1];
                    max=ac.second[i].first[0];
                }
                //std::cout <<nt <<'\t' <<min <<'\t' <<max <<std::endl;
                if(nt<min) return false;
                if(nt>max) return false;
            }
        }
        return true;
    }

    std::map<char,std::vector<int>> RotSequence::AACountRecord(
            const std::map<char,std::vector<std::pair<std::vector<int>,std::set<int>>>> &aacr){
        std::map<char,std::vector<int>> results;
        for(const auto &ac:aacr) {
            std::vector<int> ns;
            for(int i=0;i<ac.second.size();i++) {
                const std::set<int>& as = ac.second[i].second;
                int nt=0;
                for(auto & i: as) {
                    if(aaSeq[i]==ac.first) nt++;
                }
                ns.push_back(nt);
            }
            results.insert({ac.first,ns});
        }
        return results;
    }

    double DesignMC::mcRunWithAACountRestrain(RotSequence* result,
            const std::map<char,std::vector<std::pair<std::vector<int>,std::set<int>>>> &aacr){

         RotSequence unit(dt->rotGroups, aacr);
         result->copyValue(&unit);
         double energy = unit.totEnergy(dt);
         int len = unit.getLength();
         int randPos, randChoice;
         double mutEnergy;

         std::map<int,std::vector<AANumberControlInMC>> aanc;
         std::map<char,std::vector<int>> aactl=unit.AACountRecord(aacr);
         for(const auto &ac:aacr) {
             for(int i=0;i<ac.second.size();i++) {
                 const auto &p=ac.second[i];
                 AANumberControlInMC aanc1;
                 aanc1.sites = p.second;
                 aanc1.aa = ac.first;
                 aanc1.min = p.first[0];
                 aanc1.max = p.first[1];
                 aanc1.num = aactl.at(ac.first)[i];
                 for(const int &ix:p.second) {
                     if(aanc.find(ix)==aanc.end()) aanc.insert({ix,std::vector<AANumberControlInMC>()});
                     aanc.at(ix).push_back(aanc1);
                 }
             }
         }

         srand((unsigned)time(NULL));
         double minEnergy = energy;
         int count = 0;
         for(double T = T0;T>T1;T=T*annealFactor){
             bool noChange = true;
             int acNum = 0;
             for(int i=0;i<step;i++){
                 count ++;
                 randPos = rand()%len;

                 //randChoice = rand()%(unit.choiceNum(randPos));

                 char aa_old;
                 char aa_new;

                 if(aanc.find(randPos)==aanc.end()) randChoice = rand()%(unit.choiceNum(randPos));
                 else {
                     std::set<int> tobedelete;
                     char seqnow=unit.getaatype(randPos);
                     RotamerGroup* rgs=unit.getRotGroup(randPos);
                     for(AANumberControlInMC &ac:aanc.at(randPos)) {
                         if(seqnow==ac.aa) {
                             if(ac.min==ac.max || ac.min==ac.num) {
                                    for(int j=0;j<rgs->rotList.size();j++) {
                                         char c=unit.getrn().triToSin(rgs->rotList[j]->triName);
                                         if(c!=seqnow) tobedelete.insert(j);
                                    }
                             }
                         } else {
                             if(ac.min==ac.max || ac.max==ac.num) {
                                 for(int j=0;j<rgs->rotList.size();j++) {
                                     char c=unit.getrn().triToSin(rgs->rotList[j]->triName);
                                     if(c==ac.aa) tobedelete.insert(j);
                                 }
                             }
                         }
                     }

                     int nt=unit.choiceNum(randPos);
                     std::vector<int> tobechose;
                     for(int j=0;j<nt;j++) if(tobedelete.find(j)==tobedelete.end()) tobechose.push_back(j);
                     assert(tobechose.size()>0);
                     int ix = rand()%(tobechose.size());
                     randChoice = tobechose[ix];

                     aa_old = seqnow;
                     aa_new = unit.getrn().triToSin(rgs->rotList[randChoice]->triName);
                 }

                //bool isgood{false};
                 //for(int j=0;j<500;j++) {
                //    //if(unit.AACountRestrain(aacr)) {
                //         isgood=true;
                //         break;
                //     }
                //     randPos = rand()%len;
                //     randChoice = rand()%(unit.choiceNum(randPos));
                // }
                 //if(!isgood) continue;

                 mutEnergy = unit.mutEnergy(randPos,randChoice, dt);
                 if(mutEnergy < 0)
                     noChange = false;
                 if(accept(mutEnergy, T)){
                     acNum++;
                     unit.acceptMutAndUpdateSeq(randPos, randChoice);
                     energy += mutEnergy;
                     if(energy < minEnergy){
                         minEnergy = energy;
                         result->copyValue(&unit);
                     }
                     if(aa_old!=aa_new) {
                         if(aacr.find(aa_old)!=aacr.end()) {
                             for(auto &p:aacr.at(aa_old)) {
                                 auto &si=p.second;
                                 if(si.find(randPos)!=si.end()) {
                                     for(int ix:si) {
                                         for(int j=0;j<aanc.at(ix).size();j++) {
                                             if(aanc.at(ix)[j].aa==aa_old) aanc.at(ix)[j].num --;
                                         }
                                     }
                                 }
                             }
                         }
                         if(aacr.find(aa_new)!=aacr.end()) {
                             for(auto &p:aacr.at(aa_new)) {
                                 auto &si=p.second;
                                 if(si.find(randPos)!=si.end()) {
                                     for(int ix:si) {
                                         for(int j=0;j<aanc.at(ix).size();j++) {
                                             if(aanc.at(ix)[j].aa==aa_new) aanc.at(ix)[j].num ++;
                                         }
                                     }
                                 }
                             }
                         }
                     }
                     if(aanc.find(randPos)!=aanc.end()) {
                         if(aa_old!=aa_new) {
                             if(aacr.find(aa_old)!=aacr.end()) {
                                 for(int j=0;j<aacr.at(aa_old).size();j++) {

                                 }
                             }
                             if(aacr.find(aa_new)!=aacr.end()) {

                             }
                         }
                     }
                 }
             }
             printf("%7.4f %10.4f %10.4f %4d\n",T,energy,minEnergy,acNum);
             if(noChange) break;
         }
         double finalEnergy = result->totEnergy(dt);

         return finalEnergy;
    }
*/


} /* namespace NSPdesignseq */

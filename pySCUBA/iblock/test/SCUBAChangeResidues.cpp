/*
 * SCUBAChangeResidues
 *
 *  Created on: 2021326
 *      Author: wangsheng
 */
#include "iblock/molmodeler.h"
#include "iblock/intrctmol.h"
#include "iblock/intrctblck.h"
#include "designseq/StructureInfo.h"
#include "geometry/calculators.h"
#include "designseq/ProteinRep.h"
#include "dstl/randomengine.h"
#include "proteinrep/residuestate.h"
#include "proteinrep/aminoacidseq.h"
#include <dataio/splitstring.h>
#include "dataio/parameters.h"
#include <fstream>
#include <cassert>
#include <string>
#include <iostream> 

using namespace std;
using namespace NSPintrct;
using namespace NSPdesignseq;
using namespace NSPproteinrep;
using namespace NSPgeometry;

struct ChangePar {
    std::string StartPDB;
    std::string Newsequence;
    std::string Changesites;
    int MakeLVG ;
    int Changesequence ;
    int Replaceresidues ;
    int Chain = 0;
    std::string OutputFile{""};
    ChangePar() { ; }
    ChangePar(const std::vector<std::string>& controllines) {
        //for (auto s : controllines) std::cout << s << std::endl;
        std::map<std::string, double> doublepars{};
        std::map<std::string, std::vector<std::string>> stringvecpars{ };
        std::map<std::string, std::vector<double>> doublevecpars{};
        std::map<std::string, std::string> stringpars{ {"StartPDB",""},{"OutputFile",""},{"Newsequence",""}, {"Changesites","" } };
        std::map<std::string, std::vector<int>> intvecpars{};
        std::map<std::string, int> intpars{ {"MakeLVG",0},{"Changesequence",0},{"Replaceresidues",0},{"Chain",0} };
        NSPdataio::ParameterSet pset;
        pset.initdefaultkeys(doublepars, stringpars, intpars, doublevecpars, stringvecpars, intvecpars);
        pset.adjustvalues(controllines);
        pset.getval("StartPDB", &StartPDB);
        pset.getval("Newsequence", &Newsequence);
        pset.getval("MakeLVG", &MakeLVG);
        pset.getval("Changesequence", &Changesequence);
        pset.getval("Replaceresidues", &Replaceresidues);
        pset.getval("Chain", &Chain);
        pset.getval("Changesites", &Changesites);
        pset.getval("OutputFile", &OutputFile);
    }
};
/*
ChangePar::ChangePar(const std::vector<std::string>& controllines) {
    std::map<std::string, double> doublepars{};
    std::map<std::string, std::vector<std::string>> stringvecpars{};
    std::map<std::string, std::vector<double>> doublevecpars{};
    std::map<std::string, std::string> stringpars{ {"StartPDB",""},{"OutputFile",""},{"Newsequence",""} };
    std::map<std::string, std::vector<int>> intvecpars{ {"Changesites",{0,0} }};
    std::map<std::string, int> intpars{ {"MakeLVG",0},{"Changesequence",0},{"Replaceresidue",0},{"Chain",0} };
    NSPdataio::ParameterSet pset;
    pset.initdefaultkeys(doublepars, stringpars, intpars, doublevecpars, stringvecpars, intvecpars);
    pset.adjustvalues(controllines);
    pset.getval("StartPDB", &(this->StartPDB));
    pset.getval("Newsequence", &(this->Newsequence));
    pset.getval("MakeLVG", &(this->MakeLVG));
    pset.getval("Changesequence", &(this->Changesequence));
    pset.getval("Replaceresidue", &(this->Replaceresidue));
    pset.getval("Chain", &(this->Chain));
    pset.getval("OutputFile", &(this->OutputFile));
}*/


ChangePar readparameters(const std::string& parfile) {
    NSPdataio::ControlFile cf;
    cf.readfile(parfile);
    //std::cout << parfile << std::endl;
    std::string fn = "ChangePar";
    std::vector<std::string> lines = cf.getcontrolines(fn);
    //std::vector<std::string> lines = cf.getcontrolines("ChangePar");
    for (auto s : lines) std::cout << s << std::endl;
    return ChangePar(lines);
}
int main(int argc, char** argv) {
    if (argc != 2) {
        std::cout << "usage: SCUBAChangeResidues  parameterfile" << std::endl;
        exit(1);
    }

    ChangePar par = readparameters(argv[1]);
    std::cout << "read par" << std::endl;
    std::ifstream is(par.StartPDB);
    auto rstates = residuestates(is);
    int chainidx = par.Chain;
    int chaini = par.Chain;
    std::string ssseq;
    for (auto& c : rstates) {
        std::string seq;

        std::cout << "Chain " << chaini++ << ":" << std::endl;
        std::cout << "seq\tresidue\tSS\tSAI\tPHI\tPSI" << std::endl;
        int seqidx = 1;
        for (auto& r : c) {
            seq.push_back(AminoAcidSeq::name2code(r.sidechainstate.resiudetype));
            ssseq.push_back(r.backbonestate.SSState);
            std::cout << seqidx++ << "\t"
                << r.sidechainstate.resiudetype << "\t"
                << r.backbonestate.SSState << "\t"
                << r.backbonestate.sai << "\t"
                << r.backbonestate.phi << "\t"
                << r.backbonestate.psi << std::endl;
        }
        std::cout << "old residues seq: " << seq << std::endl;
        std::cout << "Secondary structure seq: " << ssseq << std::endl;
    }
    int nres = ssseq.size();






    // std::ifstream in(argv[2]);
   //  std::string seq;
   //  std::getline(in, seq);

    //std::cout << "new residues seq: " << par.Newsequence << std::endl;
    //assert(nres == par.Newsequence.size());
    std::cout << "number of residues: " << nres << std::endl;

    MolModeler modeler;
    //IntrctMol imol;
    AAConformersInModel model;
    model.readpdbfile(par.StartPDB);
    std::shared_ptr<IntrctMol> imol = std::shared_ptr<IntrctMol>(new IntrctMol(model.conformers));
    modeler.settarget(*imol);
    modeler.checkclash(true);
    if (par.Changesequence == 1) {
        std::cout << "Changesequence" << std::endl;
        std::cout << "new residues seq: " << par.Newsequence << std::endl;
        assert(nres == par.Newsequence.size());
        if (par.MakeLVG != 0 || par.Replaceresidues != 0) {
            std::cout << "cannot use multiple functions at the same time" << std::endl;
            exit(1);
        }
        for (int c = 0;c < nres;++c) {
            //if (seq[c] == 'H') {
            //    modeler.changeresidue(MolModeler::SiteIdx{ 0,c }, "LEU"); //replace residue c (starting idx=0) of chain 0 by leu
            //}
            //else {
            //    if (seq[c] == 'E') {
            //        modeler.changeresidue(MolModeler::SiteIdx{ 0,c }, "VAL");
            //    }
            //    else {
            //        modeler.changeresidue(MolModeler::SiteIdx{ 0,c }, "GLY");
            //    }
            //}
            char AA = par.Newsequence[c];
            switch (AA) {
            case'A':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "ALA");break;
            case'R':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "ARG");break;
            case'N':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "ASN");break;
            case'D':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "ASP");break;
            case'C':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "CYS");break;
            case'Q':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "GLN");break;
            case'E':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "GLU");break;
            case'G':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "GLY");break;
            case'H':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "HIS");break;
            case'I':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "ILE");break;
            case'L':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "LEU");break;
            case'K':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "LYS");break;
            case'M':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "MET");break;
            case'F':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "PHE");break;
            case'P':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "PRO");break;
            case'S':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "SER");break;
            case'T':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "THR");break;
            case'W':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "TRP");break;
            case'Y':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "TYR");break;
            case'V':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "VAL");break;

            }
        }
    }



    if (par.MakeLVG == 1) {
        std::cout << "MakeLVG" << std::endl;
        if (par.Changesequence != 0 || par.Replaceresidues != 0) {
            std::cout << "cannot use multiple functions at the same time" << std::endl;
            exit(1);
        }
        for (int c = 0;c < nres;++c) {
            if (ssseq[c] == 'H') {
                modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "LEU"); //replace residue c (starting idx=0) of chain 0 by leu
            }
            else {
                if (ssseq[c] == 'E') {
                    modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "VAL");
                }
                else {
                    modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "GLY");
                }
            }
        }
    }


    if (par.Replaceresidues == 1) {
        if (par.Changesequence != 0 || par.MakeLVG != 0) {
            std::cout << "cannot use multiple functions at the same time" << std::endl;
            exit(1);
        }
        std::cout << "replaceresidues:  " << par.Changesites << std::endl;
        std::vector<std::string> n_site;
        NSPutils::split(par.Changesites, n_site, ",");
        
        for (int n = 0;n < n_site.size();++n){
        //    std::string  onesite;

            std::vector<std::string> temp_word;
            NSPutils::split(n_site[n], temp_word, " ");
            assert(temp_word.size() == 2);
        //  onesite.push_back(temp_word[0]);

        int c = std::stoi(temp_word[0]);
        //std::cout << c << std::endl;
        std::string residue{temp_word[1]};
        char AA = residue[0];
        //modeler.changeresidue(MolModeler::SiteIdx{ 0,c }, residue);
        switch (AA) {
        case'A':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "ALA");break;
        case'R':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "ARG");break;
        case'N':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "ASN");break;
        case'D':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "ASP");break;
        case'C':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "CYS");break;
        case'Q':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "GLN");break;
        case'E':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "GLU");break;
        case'G':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "GLY");break;
        case'H':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "HIS");break;
        case'I':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "ILE");break;
        case'L':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "LEU");break;
        case'K':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "LYS");break;
        case'M':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "MET");break;
        case'F':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "PHE");break;
        case'P':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "PRO");break;
        case'S':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "SER");break;
        case'T':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "THR");break;
        case'W':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "TRP");break;
        case'Y':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "TYR");break;
        case'V':modeler.changeresidue(MolModeler::SiteIdx{ chainidx,c }, "VAL");break;

        }
        temp_word.clear();
    }
    }

    if (par.Replaceresidues != 1 && par.Changesequence != 1 && par.MakeLVG != 1) {
        std::cout << "At least one function must be specified" << std::endl;
        exit(1);
    }







    std::string filename = par.OutputFile  + ".pdb";
    std::ofstream ofs(filename);
    imol->writepdb(ofs);
    ofs.close();
    //std::string filename2 = argv[3] + std::to_string(1) + "_newseq.txt";
    //std::ofstream ofs1(filename2);
    //for (auto c : seq) ofs1 << c;
    //ofs1 << std::endl;
    //ofs1.close();


}




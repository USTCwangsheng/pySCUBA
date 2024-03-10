/*
 * @FileName:
 * @Descripttion:
 * @version:
 * @Author: Miaoyy
 * @Date: 2021-03-26 14:18:01
 * @LastEditors: Miaoyy
 * @LastEditTime: 2021-06-01 15:25:00
 */
#include "iblock/molmodeler.h"
#include "geometry/calculators.h"
#include "dstl/randomengine.h"
#include <dataio/splitstring.h>
#include <fstream>
#include <cassert>
#include <string>
#include <iostream>
#include "dataio/parameters.h"
using namespace std;
using namespace NSPintrct;
using namespace NSPproteinrep;
using namespace NSPgeometry;

struct SketchPar{
  int linkloop=1;
  std::string SketchFile;
  int randomseed=243;
  int Gennumber=1;
  std::string OutputFile;
	// std::string jobname{"task1"};
  SketchPar(){;}
	SketchPar(const std::vector<std::string> & controllines);
};
SketchPar::SketchPar(const std::vector<std::string> & controllines){
  std::map<std::string, double> doublepars {};
  std::map<std::string, std::vector<std::string>> stringvecpars {};
  std::map<std::string, std::vector<double>> doublevecpars {};
  std::map<std::string, std::string> stringpars {{"SketchFile",""},{"OutputFile",""}};
  std::map<std::string, std::vector<int>> intvecpars {};
  std::map<std::string, int> intpars{{"LinkLoop",0},{"RandomSeed",0},{"GenerationNumber",1}};
  NSPdataio::ParameterSet pset;
  pset.initdefaultkeys(doublepars,stringpars,intpars,doublevecpars,stringvecpars,intvecpars);
	pset.adjustvalues(controllines);
	pset.getval("SketchFile",&(this->SketchFile));
	pset.getval("LinkLoop",&(this->linkloop));
	pset.getval("RandomSeed",&(this->randomseed));
  pset.getval("OutputFile",&(this->OutputFile));
  pset.getval("GenerationNumber",&(this->Gennumber));
}

void readss(string file,  std::vector<char> &sseq,std::vector<std::vector<double>> &sscrd,std::vector<int>&sslen,std::vector<char> &direction){
  ifstream infile;
  // std::cout<<"start file read"<<std::endl;
  infile.open(file.data());
  // std::cout<<"file is read"<<std::endl;
  assert(infile.is_open());
  string s;
  sseq.clear();
  sscrd.clear();
  // std::cout<<"getline"<<std::endl;
  while(getline(infile,s))
  {
      std::vector<double> temp_crd_start;
      std::vector<double> temp_crd_direc;
      std::vector<int> len_range;
      std::vector<std::string> words;
      NSPutils::split(s, words, ";");
      int flag=0;
      for (auto &word : words){
        std::vector<std::string> word1s;
        NSPutils::split(word, word1s, " ");
        for (auto &word1 : word1s){
          char c=word1[0];
          if ((c == 'H' or c == 'E') && flag==0 ){
            sseq.push_back(c);
            sseq.push_back(c);
          }
          else if ((flag==0)&& (c == 'N' or c == 'C')){
            direction.push_back(c);
          }
          else if (flag==0){
            len_range.push_back(stoi(word1));
          }
          if (flag==1){
            temp_crd_start.push_back(stod(word1));
          }
          if (flag ==2){
            temp_crd_direc.push_back(stod(word1));
          }
        }
        flag++;
      }

      if (!len_range.empty())
      {
          if (len_range.size() == 2)
          {
              int len = NSPdstl::RandomEngine<>::getinstance().intrng(len_range[0], len_range[1])();
              sslen.push_back(len);
          }
          else if (len_range.size() == 1)
          {
              int len = len_range[0];
              sslen.push_back(len);
          }
          else
          {
              std::cout << "Please input according to the rules" << std::endl;
              exit(1);
          }
      }
      sscrd.push_back(temp_crd_start);
      sscrd.push_back(temp_crd_direc);

  }
  infile.close();
}


SketchPar readparameters(const std::string &parfile){
  NSPdataio::ControlFile cf;
	cf.readfile(parfile);
  std::vector<std::string> lines=cf.getcontrolines("SketchPar");
  return SketchPar(lines);
 }

int main(int argc, char**argv){
  if (argc!=2){
    std::cout <<"usage: program parameterfile" <<std::endl;
    exit(1);
  }
	SketchPar read_result= readparameters(argv[1]);

  std::vector<char> sseq;
  std::vector<std::vector<double>> sscrd;
  std::vector<int> sslen;
  std::vector<char> finalss;
  std::vector<std::string>ss_info;
  std::vector<std::string>ssEH;
  std::vector<std::string>ssE;
  std::vector<std::string>ssH;
  std::vector<std::string>ssC;
  std::vector<char> direction;
  // readss(read_result.SketchFile,sseq,sscrd,sslen,direction);
  // NSPdstl::RandomEngine<>::getinstance().reseed(read_result.randomseed);



  int i=0;
  while (i<read_result.Gennumber){
  IntrctMol imol;
  MolModeler modeler;
    finalss.clear();
    ss_info.clear();
    ssEH.clear();
    ssE.clear();
    ssH.clear();
    ssC.clear();
    readss(read_result.SketchFile,sseq,sscrd,sslen,direction);
    if (i ==0){
      NSPdstl::RandomEngine<>::getinstance().reseed(read_result.randomseed);
    }
    modeler.settarget(imol);
    modeler.buildssele(sseq,sscrd,sslen,read_result.linkloop,finalss,ss_info,  ssEH, ssE, ssH, ssC,direction);
    std::string filename = read_result.OutputFile+std::to_string(i)  + ".pdb";
    std::ofstream ofs(filename);
    imol.writepdb(ofs);
    ofs.close();
    std::string filename2=  read_result.OutputFile +std::to_string(i)+"_final_ss.txt";
    std::ofstream ofs1(filename2);
    for(auto c:finalss) ofs1<<c;
    ofs1<<std::endl;
    for(auto s:ss_info) ofs1<<s;
    ofs1 << std::endl;
    for (auto eh : ssEH) ofs1 << eh;
    ofs1 << std::endl;
    for (auto e : ssE) ofs1 << e;
    ofs1 << std::endl;
    for (auto h : ssH) ofs1 << h;
    ofs1 << std::endl;
    for (auto c : ssC) ofs1 << c;
    ofs1 << std::endl;
    ofs1.close();
    i++;
  }

}


/*
 * aasequtil.h
 *
 *  Created on: 2019年12月3日
 *      Author: hyiu
 */

#ifndef IBLOCK_AASEQUTIL_H_
#define IBLOCK_AASEQUTIL_H_
#include "designseq/ResName.h"
#include <algorithm>

/**
 * represent an amino acid sequence specified as a string of 1-letter codes
 * by a vector of strings of 3-letter code
 */
inline std::vector<std::string> getseq3l(const std::string &seq1l){
	auto &resname=NSPdesignseq::ResName::resname();
	std::vector<std::string> res;
	for(int i=0;i<seq1l.length();++i){
				res.push_back(resname.sinToTri(toupper(seq1l[i])));
	}
	return res;
}

/**
 * represent an amino acid sequence specified as a vector of 3-letter residue type
 * strings to a string of 1-letter codes
 */
inline std::string getseq1l(const std::vector<std::string> &seq3l){
	auto &resname=NSPdesignseq::ResName::resname();
	std::string res;
	for(int i=0;i<seq3l.size();++i){
			std::string r=seq3l.at(i);
			std::transform(r.begin(),r.end(),r.begin(),toupper);
			res.append(1,resname.triToSin(r));
	}
//	transform(seq3l.begin(),seq3l.end(),res.begin(),&resname.triToSin);
	return res;
}



#endif /* IBLOCK_AASEQUTIL_H_ */

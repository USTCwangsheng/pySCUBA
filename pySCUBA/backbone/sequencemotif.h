/*
 * @Descripttion: 
 * @version: 0.0
 * @Author: yuhaihan
 * @Date: 2022-03-30 10:50:45
 * @LastEditTime: 2022-04-08 10:56:20
 */
/*
 * sequencemotif.h
 *
 *  Created on: 2016年12月28日
 *      Author: hyliu
 */

#ifndef BACKBONE_SEQUENCEMOTIF_H_
#define BACKBONE_SEQUENCEMOTIF_H_

namespace NSPproteinrep{

class SequenceMotif {
public:
	template<typename ITER>
	static std::string tomotif(ITER iter,int length){
		std::string seqmotif;
		std::string omigaseq;
		for(int i=0;i<length; ++i){
			if(iter->resname=="GLY") seqmotif=seqmotif+"G";
			else if(iter->resname=="PRO"){
				seqmotif=seqmotif+"P";
			}
			else seqmotif=seqmotif+"X";
			if(iter->omiga() <-90|| iter->omiga()>90){
				omigaseq +="c";
			} else {
				omigaseq +="t";
			}
			++iter;
		}
		return seqmotif+omigaseq;
	}
};

}



#endif /* BACKBONE_SEQUENCEMOTIF_H_ */

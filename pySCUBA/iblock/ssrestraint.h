/*
 * ssrestraint.h
 *
 *  Created on: 2017年12月20日
 *      Author: hyliu
 */

#ifndef SD_SSRESTRAINT_H_
#define SD_SSRESTRAINT_H_
#include "iblock/phipsicodes.h"

namespace NSPintrct{
class SSRestraint{
public:
	SSRestraint(){;}
	SSRestraint(int chain, int sitebegin,int siteend, int ssid, double targetprob,double k):
		chainid_(chain),sitebegin_(sitebegin),siteend_(siteend),targetssid_(ssid),
		targetssratio_(targetprob),
		kres_(k){;}
	double energy(const std::vector<std::vector<SSCode>> & sscodes,
			std::vector<NSPgeometry::XYZ> *forces) const;
private:
	int chainid_{0};
	int sitebegin_{-1};
	int siteend_{-1};
	int targetssid_{-1};
	double targetssratio_{0};
	double kres_{0};
};
}



#endif /* SD_SSRESTRAINT_H_ */

/*
 * rmsdrestraintreader.h
 *
 *  Created on: 2020年9月30日
 *      Author: hyiu
 */

#ifndef RMSDRESTRAINTREADER_H_
#define RMSDRESTRAINTREADER_H_
#include <string>
#include <array>
#include <vector>
#include <memory>
namespace NSPdataio{
 class InputLines;
}
namespace  NSPintrct{
  class IntrctMol;
}
namespace NSPsampling{
class MolSDRun;
class ResTraint;
class RMSDRestraintReader{
public:
	struct RMSDTerm{
		double rmsdmin, rmsdswitch,kres;
		std::string pdbref;
		std::vector<std::pair<int,int>> refblks;
		std::vector<std::pair<int,int>> restrainedblks;
		int restrainedatoms;
};
	RMSDRestraintReader(){;};
//	void addstructrestraints(MolSDRun &sdrun);
	void readterms(NSPdataio::InputLines &input, int &lidx);
	void addrestraints(NSPintrct::IntrctMol &imol,std::vector<std::shared_ptr<ResTraint>> &results,int& switchkey);
private:
//	void addrestraint(MolSDRun  &sdrun,RMSDTerm &term);
	std::vector<RMSDTerm> terms_;
//	void readterms(const std::string &filename);
};
}
#endif /* RMSDRESTRAINTREADER_H_ */

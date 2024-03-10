/*
 * solver.h
 *
 *  Created on: 2016年12月7日
 *      Author: hyliu
 */

#ifndef LOOPCLOSURE_LOOPSOLVER_H_
#define LOOPCLOSURE_LOOPSOLVER_H_
#include "geometry/rigidtransform.h"
#include <vector>
namespace NSPloopclosure {

struct LoopSolution {
	NSPgeometry::RigidTransform rt1;
	NSPgeometry::RigidTransform rt2;
};
class LoopSolver {
public:
	enum {N1=0, CA1=1, C1=2, N2=3,CA2=4,C2=5,N3=6,CA3=7,C3=8};
	static unsigned int getSolutions(std::vector<NSPgeometry::XYZ> & atomcrds,std::vector<NSPgeometry::XYZ> &fixedcrds,
			std::vector<LoopSolution> &solutions);
	static unsigned int getSolutions(std::vector<NSPgeometry::XYZ> & atomcrds,
			std::vector<LoopSolution> &solutions);
};
}



#endif /* LOOPCLOSURE_LOOPSOLVER_H_ */

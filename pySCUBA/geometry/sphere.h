/*
 * sphere.h
 *
 *  Created on: 2020年7月23日
 *      Author: hyiu
 */

#ifndef GEOMETRY_SPHERE_H_
#define GEOMETRY_SPHERE_H_
#include "geometry/xyz.h"
namespace NSPgeometry{
struct Sphere{
		XYZ center{0,0,0};
		double radius2{1.0};
		Sphere(const XYZ &c, double r):center(c){radius2=r*r;}
		Sphere(){;}
		bool insphere(const XYZ &x) const {return ((x-center).squarednorm()<radius2);}
		template <typename RNG> XYZ randompoint(RNG &rng) const {
				 return center+XYZ(rng,sqrt(radius2));
		}
	};
}



#endif /* GEOMETRY_SPHERE_H_ */

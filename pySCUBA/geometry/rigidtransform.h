/*
 * rigidtranform.h
 *
 *  Created on: 2016年11月8日
 *      Author: hyliu
 */

#ifndef RIGIDTRANSFORM_H_
#define RIGIDTRANSFORM_H_

#include "geometry/rotation.h"

namespace NSPgeometry {
/**A rigidtransform comprises a rotation and a translation.
 *
 * X_transformed=rotation(X_original)+translation
 */
class RigidTransform{

public:

	RigidTransform(){;}
	/** constructor
	 *
	 * @param Q: the roation represented as a quternion
	 * @param rotcenter : a point on the rotation axis
	 * @param t the translation
	 */
	RigidTransform(const QuaternionCrd &Q, const XYZ &rotcenter=XYZ(), const XYZ &t=XYZ()):trans_(t){
		rot_.init(Q,rotcenter);
	}
	/**
	 *reformulate a rotation into  a rotation around the origin and a translation
	 */
	RigidTransform(const Rotation &a){
		for(int i=0;i<3;i++) {
			for(int j=0;j<3;j++ )
			rot_.matrix()[i][j]=a.matrix()[i][j];
		}
		trans_=a.center()-rot_.applytoCopy(a.center());
	}

	/** constructor
	 *
	 * @param r: the rotation
	 * @param t: the translation
	 */
	RigidTransform(const Rotation & r, const XYZ & t): rot_(r),trans_(t){;}

	RigidTransform axistoorigin() const {
		RigidTransform rt0(rot_);
		rt0.trans_ = rt0.trans_+trans_;
		return rt0;
	}

	/**
	 * The constructed rigidtransform, when applied to rl, will set the internal coordinates of point
	 * rl with respect to point ri, rj and rk to given values.
	 * ri, rj, rk must not be colinear.
	 */
	RigidTransform(const XYZ & ri, const XYZ &rj, const XYZ &rk, const XYZ &rl,
			double dkl, double tjkl, double pijkl);

	RigidTransform getreverse() const {
		RigidTransform r0=this->axistoorigin();
		RigidTransform res;
		for(int i=0;i<3;++i)
			for(int j=0;j<3;++j)	res.rotation().matrix()[i][j]=r0.rotation().matrix()[j][i];
		res.translation()=-res.applytoCopy(r0.translation());
		return res;
	}
	/**
	 * Obtain a rigidtransform to be applied to a set of moving points so that their
	 * positions relative to a set of fixed points are changed in defined ways.
	 * A joint is defined by three points from the fixed set, ifx,jfx,and kfx, together
	 * with the moving set, kmv,jmv,and imv.  The  Cartesian coordinates of these points
	 * are given in jointpts in the order of (ifx,jfx,kfx,kmv,jmv,imv).
	 * The target relative geometries are given by the following internal coordinate values
	 * in intcrds.
	 * torsion(ifx,jfx,kfx,kmv),angle(jfx,kfx,kmv),
	 * distance(kfx,kmv),torsion(jfx,kfx,kmv), angle(kfx,kmv,jmv) and torsion(kfx,kmv,jmv,imv).
	 *Applying the constructed Rigidtransform to kmv,jmv and imv will lead to these given internal
	 *coordinate values.
	 */
	RigidTransform(const std::vector<XYZ> &jointpts, const std::vector<double> &intcrds);

	/** apply the rigid transform to a point
	 *
	 * @param p: point to the XYZ point to be transformed
	 */

	void apply(XYZ *p) const;

	/**Apply the rigid transform to a copy of an input point
	 *
	 * @param p: input XYZ point, will not be changed.
	 * @return the transformed XYZ point.
	 */
	XYZ applytoCopy(const XYZ & p) const;
	/**get the rotation part of the rigid transform
	 *
	 */
	Rotation &rotation() {return rot_;}
	/** get the const rotation part of the rigid transform
	 *
	 */
	const Rotation & rotation() const {return rot_;}
	/** get the translation part of the rigid transform
	 *
	 */
	XYZ & translation() {return trans_;}
	/** get the const translation oart of the rigid transform
	 *
	 */
	const XYZ & translation() const {return trans_;}
private:
	/**the rotation
	 *
	 */
	Rotation rot_;

	/**The translation
	 *
	 */
	XYZ trans_;
};
/**Combine a rigid transform with a rotation, get a new rigid transform
 * X_transformed=Rotation_r(RigidTranform rt(X_original)
 * @return the resulting rigid transform
 */
RigidTransform applyRotation(const Rotation & r, const RigidTransform &rt);
/**Generate a rigid transformation with random roation and translation
 * The rotation is around a random axis passing the origin, with uniform distributed angles between
 * 0 and maxrotate(in degrees)
 * The translation is of length between 0 to maxtranslate
 *
 */
RigidTransform randomrigidtransform(double maxrotate, double maxtranslate);

/**Get the reigidtransformation that would superimpose the coordinate set crdb onto the
 * coordinate set crda.
 * @param crda: the reference(fixed) coordinate
 * @param crdb: the coordinates to be superimposed (only the transform is calculated, the
 * coorinates are not changed.
 * @param alignedpositions: the aligned positions between crda and crdb. Only the atoms in the aligned positions will be
 * considered for superposition.
 * @param dev2: if not null, the squared rmsd will be delivered to dev2.
 *
 */
RigidTransform superpose(const std::vector<XYZ> & crda, const std::vector<XYZ> &crdb,
		const std::vector<std::pair<int,int>> & alignedpostions,double *dev2=nullptr);

RigidTransform operator *(const RigidTransform &  a, const RigidTransform &b);
inline RigidTransform operator *(const Rotation &a, const Rotation &b){
		return RigidTransform(a) *RigidTransform(b);
}
}
#endif /* RIGIDTRANSFORM_H_ */

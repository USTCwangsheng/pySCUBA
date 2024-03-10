/*
 * rectgrid.h
 *
 *  Created on: 2020年7月22日
 *      Author: hyiu
 */

#ifndef GEOMETRY_RECTGRID_H_
#define GEOMETRY_RECTGRID_H_
#include "geometry/xyz.h"
#include <map>
namespace NSPgeometry{
struct RectGrid{
	XYZ p000; //<*left-front-bottom point
	double step;
	std::array<int,3> sizes;
	typedef std::array<int,3> PntIdx;
	RectGrid(const XYZ &center, const std::array<int,3> & gridsizes, double gridstep):
	  sizes(gridsizes),step(gridstep){
		for(int i=0;i<3;++i){
			p000[i]=center[i]-step*0.5*((double)sizes[i]-1.0);
		}
	}
	RectGrid(const XYZ &p0,const XYZ &p1,double gridstep):p000(p0),step(gridstep){
	      for(int i=0;i<3;++i) sizes[i]=(p1[i]-p0[i]-0.001*step)/step+1;
	}
	XYZ center() const{
		XYZ c;
		for(int i=0;i<3;++i) c[i]=p000[i]+step*0.5*((double) sizes[i]-1.0);
		return c;
	}
	XYZ idx2xyz(const PntIdx & pidx) const{
		XYZ p;
		for(int i=0;i<3;++i) p[i]=p000[i]+step*(double) pidx[i];
		return p;
	}
	double mindis2(const XYZ &p, const PntIdx &cidx) const{
		XYZ diff;
		XYZ c=idx2xyz(cidx);
		double hstep=0.5*step;
		for(int i=0;i<3;++i){
			double d=p[i]-c[i];
			if(d<-hstep) diff[i]=-d+hstep;
			else if(d>hstep) diff[i]=d-hstep;
			else diff[i]=0.0;
		}
		return diff.squarednorm();
	}
	std::vector<PntIdx> neighborgridpoints(XYZ c, double rcut2) const{
		 PntIdx cp,pmin,pmax;
		 int buff=rcut2/(step*step)+1;
		 for(int i=0;i<3;++i) {
			 cp[i]=(c[i]-p000[i])/step;
			 pmin[i]=cp[i]-buff;
			 pmax[i]=cp[i]+buff;
		 }
		 std::vector<PntIdx> nbs;
		 PntIdx p;
		 for(int p0=pmin[0];p0<pmax[0]+1;++p0){
				 if(p0<0 || p0>=sizes[0]) continue;
				 p[0]=p0;
				 for(int p1=pmin[0];p1<pmax[0]+1;++p1){
					  if(p1<0 || p1>=sizes[1]) continue;
					  p[1]=p1;
					  for(int p2=pmin[0];p2<pmax[0]+1;++p2){
						  if(p2<0 || p2>=sizes[2]) continue;
						  p[2]=p2;
						  if (p0==pmin[0] || p0==pmax[0] ||
								  p1==pmin[1] ||p1==pmax[1] ||
								  p2==pmin[2] ||p2==pmax[2]){
							   if(mindis2(c,p)>rcut2) continue;
						  }
						  nbs.push_back(p);
					  }
				 }
			 }
		 return nbs;
		 }
	std::map<PntIdx,std::vector<int>> neighborxyzs(const std::vector<XYZ> &crds,
			double rcut2){
		std::map<PntIdx,std::vector<int>> res;
		bool manypnts=crds.size()>1000;
		if(manypnts){
			PntIdx pidx;
			for(int i=0;i<sizes[0];++i){
				pidx[0]=i;
				for(int j=0;j<sizes[1];++j){
					pidx[1]=j;
					for(int k=0;k<sizes[2];++k){
					   pidx[2]=k;
					   res[pidx]=std::vector<int>();
					}
				}
			}
		}
		for(int i=0;i<crds.size();++i){
			auto nbs=neighborgridpoints(crds[i],rcut2);
			for(auto &pidx:nbs){
				if(manypnts) {
					res[pidx].push_back(i);
					continue;
				}
				auto entry=res.find(pidx);
				if(entry== res.end()){
					res[pidx]=std::vector<int>();
					res[pidx].push_back(i);
				} else{
					entry->second.push_back(i);
				}
			}
		}
		return res;
	}
};

}

#endif /* GEOMETRY_RECTGRID_H_ */

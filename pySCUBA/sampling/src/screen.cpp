/*
 * screen.cpp
 *
 *  Created on: 2021年6月4日
 *      Author: yxchen
 */

#include "iblock/intrctmol.h"
#include "iblock/intrctblck.h"
#include "sampling/screen.h"
#include <fstream>
#include <sstream>
#include <dirent.h>
#include "proteinrep/aaconformer.h"

using namespace std;
using namespace NSPintrct;
using namespace NSPproteinrep;
using namespace Screen;

/*
 * used for SCUBAHBScreen and SCUBACavityScreen.
 */

bool Screen::expose(XYZ crd, string cid, int rid, vector<vector<AAConformer>> conformers, double expose_th, string obj2)
{
	// collect.
	vector<XYZ> insph; // atms within 4.5A sphere.
	for (auto c : conformers)
		for (auto r : c)
			for (auto a : r.atomlist)
			{
				if (a[0] == 'H') continue;
				if (obj2 == "MC")
				{
					if (a != "N" && a != "CA" && a != "C" && a != "O") continue;
				}
				if (obj2 == "SC")
				{
					if (a == "N" || a == "CA" || a == "C" || a == "O") continue;
				}
				XYZ co = r.globalcrd[a];
				double dis = co.distance(crd);
				if (dis <= 5.4 && dis > 0)
					insph.push_back(co);
			}
	set<int> tp; // taken points
	SpherePoints sp; // local crd for 2.7A sphere points.
	for (int p = 0; p < sp.points.size(); p++)
	{
		bool cover = false;
		for (auto is : insph)
		{
			if (is.distance(crd+sp.points[p]) <= 2.7) cover = true;
			if (cover)
			{
				tp.insert(p);
				break;
			}
		}
	}
	/*
	// projection.
	set<int> tp; // taken points by projection;
	SpherePoints sp; // local crd for 4.5A sphere points.
	for (auto is : insph)
		for (int p = 0; p < sp.points.size(); p++)
		{
			if (tp.count(p) > 0) continue;
			auto pcrd = sp.points[p];
			XYZ ci = is-crd;
			double coscpci = ci*pcrd/(4.5*ci.length());
			double dis = ci.length()*sin(acos(coscpci));
			if (dis <= 0.7) tp.insert(p);
 		}
 	*/
	ofstream ofsd("expose_detailed.txt", ios::app);
	ofsd << " expose: " << cid << " " << rid << " " << 1-double(tp.size())/sp.points.size() << endl;
	if (double(tp.size())/sp.points.size() < 1-expose_th) return true;
	else return false;
}
bool Screen::formhb(XYZ crd, string cid, int rid, vector<vector<AAConformer>> conformers, double hb_th, string obj2)
{
	for (auto c : conformers)
		for (auto r : c)
		{
			if (r.chainid_or == cid[0] && abs(r.residueid_or-rid) <= 1) continue; // itself and neighbor donot form hb.
			for (auto it = r.globalcrd.begin(); it != r.globalcrd.end(); it++)
			{
				if (obj2 == "MC")
				{
					if (it->first != "N" && it->first != "O") continue;
				}
				if (obj2 == "SC")
				{
					if (it->first == "N" || it->first == "O") continue;
				}
				if (it->first[0] == 'N' || it->first[0] == 'O')
					if (crd.distance(it->second) < hb_th)
						return true;
			}
		}
	return false;
}

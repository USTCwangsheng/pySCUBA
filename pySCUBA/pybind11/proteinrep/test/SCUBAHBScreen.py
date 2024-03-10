#!/bin/python3
# This program calculate the score of MC hygrogen bonds for a protein
# ori_cpp: /sampling/test/SCUBAHBSCreen.cpp
import sys,os
path='/home/cyx/workspace/pySCUBA/pybind11/bin/'
sys.path.append(path)
import pb_geometry
import pb_proteinrep
import pb_screen
"""This program calculate the score of MC hygrogen bonds for a protein"""

# input parameters:
pdb = 'start.pdb'
expose_th = 0.14
HB_dist = 3.5

# program calculations
nexp_patm = 0.0
nexp_nhbatm = 0.0
AAC = pb_proteinrep.AAConformersInModel()
AAC.readpdbfile(pdb)
cons = AAC.conformers
for c in cons:
    for r in c:
        crds = r.globalcrd
        if not(pb_screen.SpherePoints.expose(r.globalcrd["N"], r.chainid_or, r.residueid_or, cons, expose_th, 'MC')):
            nexp_patm+=1
        if not(pb_screen.SpherePoints.formhb(r.globalcrd["N"], r.chainid_or, r.residueid_or, cons, HB_dist, 'MC')):
            nexp_nhbatm+=1
        if not(pb_screen.SpherePoints.expose(r.globalcrd["O"], r.chainid_or, r.residueid_or, cons, expose_th, 'MC')):                                                       nexp_patm+=1
        if not(pb_screen.SpherePoints.formhb(r.globalcrd["O"], r.chainid_or, r.residueid_or, cons, HB_dist, 'MC')):
            nexp_nhbatm+=1
nhb = nexp_nhbatm/nexp_patm
print("score: %f" %nhb)
# score for 40 native scaffolds is 0.08957

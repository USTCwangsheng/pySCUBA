#!/bin/python3
# This program calculate rmsd between two pdbs for their ALL/MC/CA atoms.
# ori_cpp: no original cpp program.
import sys,os
path='/home/cyx/workspace/pySCUBA/pybind11/bin/'
sys.path.append(path)
import pb_geometry
import pb_proteinrep
"""This program calculate rmsd between two pdbs for their ALL/MC/CA atoms."""
# input parameters:
pdb1 = 'output.pdb'
pdb2 = 'start.pdb'
rmsdtype = 'ALL'
# rmsdtype could be ALL MC CA

# program calculations
AAC1 = pb_proteinrep.AAConformersInModel()
AAC1.readpdbfile(pdb1)
AAC2 = pb_proteinrep.AAConformersInModel()
AAC2.readpdbfile(pdb2)

def getcrd_by_name(aac,atmname=None):
    res_lst=[r.globalcrd for c in aac.conformers for r in c]
    if atmname:
        return [res[i] for res in res_lst for i in res.keys() if i in atmname]
    else:
        return [res[i] for res in res_lst for i in res.keys()]

rmsdtype_lst={"ALL":None,"MC":("N","CA","C","O"),"CA":("CA")}
atmname=rmsdtype_lst[rmsdtype]
crd1,crd2=getcrd_by_name(AAC1),getcrd_by_name(AAC2)
rmsd=pb_geometry.QuatFit.rmsd(crd1,crd2)
print("rmsd: %f" %rmsd)

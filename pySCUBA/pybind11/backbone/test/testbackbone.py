import sys
import os
# sys.path.append("/home/zhanglu/workspace/new/pybind11/bin")
sys.path.append("/home/zhanglu/workspace/pyscuba/pySCUBA/pybind11/bin")
import pb_backbone

obj_backbone=pb_backbone.BackBoneSite
# obj_backbone.resname="ALA"
# print(obj_backbone.resname)


obj_backbonebuilder=pb_backbone.BackBoneBuilder

#readbackbonefrompdb
# onebackbone=pb_backbone.readbackbonefrompdb('/home/zhanglu/data/comA_rapF/new/3/1/3ulq_papersites_ori_longhelix.pdb')
# print('len: '+str(len(onebackbone)))
# for n in range(len(onebackbone)):
#     for m in range(len(onebackbone[n])):
#         print(onebackbone[n][m].resname)

#test generaterandombackbone writeSitesToPDB
helixregion=[]
strandregion=[]
newbb=pb_backbone.generaterandombackbone(10,helixregion,strandregion)
file='/home/zhanglu/workspace/pyscuba/pySCUBA/pybind11/backbone/test/test.pdb'
pb_backbone.writeSitesToPDB(file,newbb)
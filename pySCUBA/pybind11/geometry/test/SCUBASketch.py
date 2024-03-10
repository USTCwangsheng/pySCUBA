# iblock/test/SCUBASketch.cpp
import sys
sys.path.append("/home/yuhaihan/workspace/test_py_SCUBA/pySCUBA/pybind11/bin/")
from pb_geometry import *
from pb_iblock import *
from pb_dstl import *

read_result=SketchPar.readparameters("sketch.par")

RandomEngine.reseed(read_result.randomseed)

for i in range(read_result.Gennumber):
    imol=USRIntrctMol()
    modeler=MolModeler()
    mainpar=SCUBASketchMainPAR()
    mainpar.readss(read_result.SketchFile)
    modeler.settarget(imol)
    modeler.buildssele(mainpar.sseq,mainpar.sscrd,mainpar.sslen,read_result.linkloop,mainpar.finalss,mainpar.ss_info,mainpar.ssEH,mainpar.ssE,mainpar.ssH,mainpar.ssC,mainpar.direction)
    imol.usrwritepdb("ouput_{}.pdb".format(i))


from platform import libc_ver
import sys
lib_path="/home/yuhaihan/workspace/test_py_SCUBA/pySCUBA/pybind11/bin/"
sys.path.append(lib_path)
from pb_geometry import *
from pb_iblock import *
from pb_dstl import *
import numpy as np

class Sketch:
    def __init__(self,sketchfile=None) -> None:
        self.random_linkloop=False
        self.para_txt=[]
        if sketchfile:
            self.read_parfile(sketchfile)
    def __repr__(self):
        return "".join(self.para_txt)
    def draw_3D(self):
        pass
    def read_parfile(self,sketchfile):
        with open(sketchfile,'r') as f:
            self.para_txt=f.readlines()
        self.mainpar=SCUBASketchMainPAR()
        self.mainpar.readss(sketchfile)
    def add_element(self,sstype,int_range,start_from_N_or_C,start_xyz,direction):
        assert int_range[0]<=int_range[1]
        assert len(start_xyz)==3
        assert len(direction)==3
        direction=np.array(direction)/np.linalg.norm(direction)
        str_="{} {} {} {};{} {} {}; {} {} {}\n".format(sstype,int_range[0],int_range[1],start_from_N_or_C,*start_xyz,*direction)
        self.para_txt.append(str_)
    def save_parfile(self,outputfile="saved_SCUBASketch.par"):
        print("the file {} will be create or replaced".format(outputfile))
        with open(outputfile,"w") as f:
            f.writelines(self.para_txt)
    def build(self):
        mainpar=self.mainpar
        return (mainpar.sseq,mainpar.sscrd,mainpar.sslen,self._link_loopmethod(),mainpar.finalss,mainpar.ss_info,mainpar.ssEH,mainpar.ssE,mainpar.ssH,mainpar.ssC,mainpar.direction)
    def pop(self):
        if self.para_txt:
            self.para_txt.pop()
        else:
            print("The sketch is clean")
    def _link_loopmethod(self):
        return 0 if self.random_linkloop else 1
    
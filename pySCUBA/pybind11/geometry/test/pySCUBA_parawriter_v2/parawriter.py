# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 09:43:11 2022

@author: yu
"""

from collections import defaultdict
# from tempfile import NamedTemporaryFile
import pickle
import copy
from usr_infowriter import *

def strcat(str_):
    return "".join(["\t=\t",str_,"\n"])
def simple(value):
    return strcat(str(value))
def values(v_lst):
    str_lst=[str(v) for v in v_lst]
    return "{}\t"[:-1].format(*str_lst)
def process(value):
    if isinstance(value,bool):
        value=1 if value else 0
    return str(value)

key_f=defaultdict(lambda:simple,{"AnnealingScheme_":values,"StoreTopConfig_":values})

class Parawriter:
    def __init__(self,fname,objs) -> None:
        self.fname=fname
        self.blocks=[]
        for obj in objs:
            self.blocks.append(Block(getattr(objs,obj)))

class Block:
    def __init__(self,obj):
        self.name=obj.__class__.__name__
        self.context=["START {}\n".format(self.name)]
        for attr in obj:
            self.add(attr,getattr(obj,attr))
        self.finish()
    def add(self,key,value):
        key=str(key)
        self.context.append(key+key_f[key](process(value)))
    def finish(self):
        self.context.append("END {}\n".format(self.name))
        
class USR_interface:
    def __init__(self,default_file="default_para.pickle",loadfiles=["usr_para.pickle"]):
        self.paras={}
        self._default_paras = {}
        self.load_default(default_file)
        if loadfiles:
            self.load_files(loadfiles)
    def load_default(self,file):
        with open(file,"rb") as f:
            py_data=pickle.load(f)
        for strategy in py_data:
            self._default_paras[strategy]=py_data[strategy]
            self.paras[strategy]=py_data[strategy]
    def load_file(self,file):
        try:
            with open(file,"rb") as f:
                py_data=pickle.load(f)
        except :
            print("file {} is empty! OR other error".format(file))
            return 0
        for strategy in py_data:
            if strategy in self.paras:
                print("Warning !Can not reads paras named {} in file {}\
                      because same name exits".format(strategy,file))
                continue
            self.paras[strategy]=py_data[strategy]
    def load_files(self,loadfiles):
        for file in loadfiles:
            self.load_file(file)
    def rewrite_strategy(self,strategy_lst,fname):
        with open(fname,'wb') as f:
            self._edit_file(strategy_lst,f)
    def add_strategy(self,strategy_lst,fname):
        pass
    def _edit_file(self,strategy_lst,f):
        for sname in strategy_lst:
            self.name_exits(sname)
        py_data={sname:value for sname,value in self.paras.items() if sname in strategy_lst}
        pickle.dump(py_data,f)
    def show_all(self):
        count=1
        for s in self.paras:
            print("{}.\t{}\n\tDescription:\t{}".format(count,s,self.paras[s].description))
            count+=1
    def copy_strategy(self,old_strategy_name,new_strategy_name):
        if self.name_exits(old_strategy_name):
            if not self.name_exits(new_strategy_name,ck_old=False):
                self.paras[new_strategy_name]=copy.deepcopy(self.paras[old_strategy_name])
    def name_exits(self,strategy_name,ck_old=True) ->bool :
        t=True if strategy_name in self.paras else False
        if ck_old and not t:
            print("name Error : strategy named {} not exits".format(strategy_name))
        if not ck_old and t:
            print("name invalid: strategy named {} already exits".format(strategy_name))
        return t
    def edit_strategy(self,strategy):
        if self.name_exits(strategy):
            return self.paras[strategy]
    def write_parafile(self,strategy,fname):
        writer=Parawriter(fname,self.paras[strategy])
        with open(fname,'w') as f:
            for block in writer.blocks:
                f.writelines(block.context)
if __name__=="__main__":
    t=USR_interface()
    t.copy_strategy("default_strategy","usr_strategy")
    t.show_all()
    usr_para=t.edit_strategy("usr_strategy")
    usr_para.SDRunPar.SavePDBSteps=100
    t.rewrite_strategy(["usr_strategy"],"usr_para.pickle")
    t.write_parafile("usr_strategy", 'sd.par')

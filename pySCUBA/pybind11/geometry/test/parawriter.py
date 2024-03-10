# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 09:43:11 2022

@author: yu
"""

from collections import defaultdict


def strcat(str):
    return "".join(["\t=\t",str,"\n"])
def simple(value):
    return strcat(str(value))
def values(v_lst):
    str_lst=[str(v) for v in v_lst]
    return "{}\t"[:-1].format(*str_lst)

key_f=defaultdict(lambda:simple,{"AnnealingScheme":values,"StoreTopConfig":values})

class Parawriter:
    def __init__(self,filename) -> None:
        self.fname=filename
        self.blocks={}
        pass
    def start_end(self,block_name):
        for block in self.blocks:
            self.blocks[block].context.append()

class Block:
    def __init__(self,blk_name) -> None:
        self.name=blk_name
        self.context=["START {}\n".format(self.name)]
    def add(self,key,value):
        key=str(key)
        self.context.append(key+key_f[key])
    def finish(self):
        self.context.append("END {}\n".format(self.name))
        



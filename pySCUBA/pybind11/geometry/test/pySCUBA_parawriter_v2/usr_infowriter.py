# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 20:46:24 2022

@author: yu
"""
import pickle
from copy import deepcopy

class Structure:
    # Class variable that specifies expected fields
    _fields= []
    def __init__(self, *args, **kwargs):
        if len(args) != len(self._fields):
            raise TypeError('Expected {} arguments'.format(len(self._fields)))
        # Set the arguments
        for name, value in zip(self._fields, args):
            setattr(self, name, value)
        # Set the additional arguments (if any)
        extra_args = kwargs.keys() - self._fields
        for name in extra_args:
            setattr(self, name, kwargs.pop(name))
        if kwargs:
            raise TypeError('Duplicate values for {}'.format(','.join(kwargs)))
    #return a iter of stirng of class keys
    def __iter__(self):
        fields=deepcopy(self._fields)
        return iter(fields)
    def keys(self):
        return str(self._fields)

class checkedmeta(type):
    def __new__(cls, clsname, bases, methods):
        # Attach attribute names to the descriptors
        for key, value in methods.items():
            if isinstance(value, Descriptor):
                value.name = key
        return type.__new__(cls, clsname, bases, methods)


# Base class. Uses a descriptor to set a value
class Descriptor:
    def __init__(self, name=None, **opts):
        self.name = name
        for key, value in opts.items():
            setattr(self, key, value)
    def __set__(self, instance, value):
        instance.__dict__[self.name] = value
        
# Descriptor for enforcing types
class Typed(Descriptor):
    expected_type = type(None)
    def __set__(self, instance, value):
        if not isinstance(value, self.expected_type):
            raise TypeError('expected ' + str(self.expected_type))
        super().__set__(instance, value)
        
# Descriptor for enforcing values
class Unsigned(Descriptor):
    def __set__(self, instance, value):
        if value < 0:
            raise ValueError('Expected >= 0')
        super().__set__(instance, value)

# Descriptor for given length
class MaxSized(Descriptor):
    def __init__(self, name=None, **opts):
        if 'size' not in opts:
            raise TypeError('missing size option')
        super().__init__(name, **opts)
    def __set__(self, instance, value):
        if len(value) >= self.size:
            raise ValueError('size must be < ' + str(self.size))
        super().__set__(instance, value)

class Integer(Typed):
    expected_type = int
class UnsignedInteger(Integer, Unsigned):
    pass
class Float(Typed):
    expected_type = float
class UnsignedFloat(Float, Unsigned):
    pass
class String(Typed):
    expected_type = str
class SizedString(String, MaxSized):
    pass
class Bool(Typed):
    expected_type = bool


class MolSystmPar(Structure,metaclass=checkedmeta):
    ChangePDBSequence	=	Bool()
    PrintParameters	=	Bool()
    JobName	=	String()
    PDBStart	=	String()
    Sequence1L	=	String()
    _fields=["ChangePDBSequence","PrintParameters","ActiveResidues",
            "FixedResidues","JobName","MainChainFixedResidues",
            "PDBStart","SideChainActiveResidues","SoftSideChainResidues",
            "Sequence1L"]
default_MolSystmPar=MolSystmPar(False,True,'',
                                '','sd0','',
                                '','',''
                                ,'')

class InteractionPar(Structure,metaclass=checkedmeta):
    PrintParameters	=	Bool()
    WriteEneDetails	=	Bool()
    CovalentWeight	=	UnsignedFloat()
    LocalHBWeight	=	UnsignedFloat()
    LocalStrWeight	=	UnsignedFloat()
    MCStericWeight	=	UnsignedFloat()
    PhiPsiWeight	=	UnsignedFloat()
    RotamerWeight	=	UnsignedFloat()
    SCPackingWeight	=	UnsignedFloat()
    SitePairWeight	=	UnsignedFloat()
    JobName	        =   String()
    _fields=["PrintParameters","WriteEneDetails","CovalentWeight",
            "LocalHBWeight","LocalStrWeight","MCStericWeight",
            "PhiPsiWeight","RotamerWeight","SCPackingWeight",
            "SitePairWeight","JobName"]
default_InteractionPar=InteractionPar(True,False,1.0,
                        0.6,0.5,1.0,
                        2.0,2.4,3.1,
                        0.32,'')


class SDRunPar(Structure,metaclass=checkedmeta):
    DOAnnealing	=	Bool()
    DOShake	=	Bool()
    PrintParameters	=	Bool()
    PrintSteps	=	UnsignedInteger()
    RandomSeed	=	UnsignedInteger()
    RecalcNeighborListSteps	=	UnsignedInteger()
    RecalcSSSteps	=	UnsignedInteger()
    SavePDBSteps	=	UnsignedInteger()
    GAMMA	=	UnsignedFloat()
    RestraintsFile	=	String()
    TimeStep	=	UnsignedFloat()
    GroupTemperatures	=	UnsignedFloat()
    AnnealingGroup	=	Integer()
    JobName	=	String()
    OutPDBFile	=	String()
    _fields=["DOAnnealing","DOShake","PrintParameters",
            "StoreTopConfig","PrintSteps","RandomSeed",
            "RecalcNeighborListSteps","RecalcSSSteps","SavePDBSteps",
            "GAMMA","RestraintsFile","TimeStep",	     
            "AnnealingScheme","GroupTemperatures","AnnealingGroup",
            "JobName","OutPDBFile","TemperatureGroups"]
default_SDRunPar=SDRunPar(True,True,True,
                          "50 2 0.5 0",50,36,
                          50,500,100,
                          1.0,"",2e-3,
                          '2.5    0.5    10000 4000 1000',0.5,-1,
                          "","finall_out1.pdb","ALL")

class Strategy(Structure):
    _fields=["MolSystmPar","InteractionPar","SDRunPar"]
    
default_strategy=Strategy(default_MolSystmPar,default_InteractionPar,default_SDRunPar, 
  description="This strategy is used for Global Stochastic Dynamic Simulated Annealing Sampling of Protein Structure")

if __name__=="__main__":
    paras={"default_strategy":default_strategy}
    with open("default_para.pickle",'wb') as f:
        pickle.dump(paras, f)
    
    
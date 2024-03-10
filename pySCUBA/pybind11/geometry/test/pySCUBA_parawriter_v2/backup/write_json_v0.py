# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 14:44:43 2022

@author: yu
"""
import json

default_para={
"default_strategy":{
    "_description":"This strategy is used for .....",
    "MolSystmPar":{
        "ChangePDBSequence"	:	0,
        "PrintParameters"	:	0,
        "ActiveResidues"	:	"",
        "FixedResidues"	:"chain0 0, 7-15 , 28-32",
        "JobName"	:   "sd0",
        "MainChainFixedResidues"	:	"",
        "PDBStart"	:  "start.pdb",
        "SideChainActiveResidues"	:	"",
        "SoftSideChainResidues"	:	"",
        "Sequence1L"	:""},
    "InteractionPar":{
        "PrintParameters"	:	0,
        "WriteEneDetails"	:	0,
        "CovalentWeight"	:	1.0,
        "LocalHBWeight"	:	0.6,
        "LocalStrWeight"	:	0.5,
        "MCStericWeight"	:	0.01,
        "PhiPsiWeight"	:	2.0,
        "RotamerWeight"	:	2.4,
        "SCPackingWeight"	:	3.1,
        "SitePairWeight"	:	0.5,
        "JobName"	:""	},
    "SDRunPar":{
        "DOAnnealing"	:	1,
        "DOShake"	:	1,
        "PrintParameters"	:	1,
        "StoreTopConfig"	:	"50	2	0.5	0",
        "PrintSteps"	:	50,
        "RandomSeed"	:	36,
        "RecalcNeighborListSteps"	:	50,
        "RecalcSSSteps"	:	500,
        "SavePDBSteps"	:	100,
        "GAMMA"	:	0.5,
        "RestraintsFile" : "restraints.txt",
        "TimeStep"	:	0.002,
        "AnnealingScheme"	:	"3    0.5    10000 4000 1000",
        "GroupTemperatures"	:	1.0,
        "AnnealingGroup" : -1,
        "JobName"	:	"",
        "OutPDBFile"	:	"finall_out1.pdb",
        "TemperatureGroups"	:	"ALL",}
    }
}
# json_str = json.dumps(default_para)
# with open ("default_para.json","w") as f:
#     json.dump(json_str,f)

with open ("default_para.json","r") as f:
    js_data=json.load(f)
    
py_data=json.loads(js_data)
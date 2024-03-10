import math
import random
import sys,os
path='/home/wangsheng/pySCUBA/pybind11/bin/'
sys.path.append(path)
from pb_iblock import *
from pb_sampling import *
from newclassscubalib import *
from pb_designseq import *

kbt=2.4942
A2NM=0.1
#################sketch生成初始结构#####################
sketchpar='newsketchpar.txt'#参数文件，没有会自动生成
skecthfile='3pvh.txt'       #sketch坐标文件 。必须
linkoption=1                #是否连接loop，0不连，1连，默认0
randomseed=15               #随机数种子
outputname='build'          #输出pdb的文件名,以此为前缀：build0.pdb,并会生成SSinfo：build0_final_ss.txt  。必须
gennumber=2                 #生成的结构数目，默认1

Sketchbuildpdb(sketchpar,skecthfile,linkoption,randomseed,outputname,gennumber)


#############SD##############
sdpdb='build0.pdb'
SSinfo='build0_final_ss.txt'
restraints='restraints1.txt'
# writerestraints(sdpdb,SSinfo,restraints)  #只自动生成rmsd，helix，scale，rg restraints，不包含contact restraints，推荐自己编辑
mpara,ipara,spara=SCUBASDinit(sdpdb,restraints)  # mpara: Molecular parameters, ipara: interaction parameter, spara: Molecular dynamics parameters
#可修改的参数展示
mpara.pdbstart = sdpdb
mpara.jobname = 'sd0'
ipara.enedetails = 0
ipara.weight_steric = 0.01
ipara.weight_phipsi = 2.0
ipara.weight_localstr = 0.5
ipara.weight_sitepair = 0.32
ipara.weight_rotamer = 2.4
ipara.weight_scpacking = 3.1
ipara.weight_localhb = 0.6
spara.annealinggroup = -1
spara.annealingscheme = [2, 0.5, 10000, 4000, 1000]
spara.atomgroups = []
spara.atommass = 14.0
spara.doannealing = 1
spara.doshake = 1
spara.emaxforstore = -1000
spara.epot_decay = 0.5
spara.jobname = ''
spara.newnbliststeps = 50
spara.newsscodesteps = 500
spara.nr_rmsd = 0.1   #单位纳米*埃,1nm=10A
spara.outpdbfile = 'sdout.pdb'
spara.topnpdbname = 'enetop_'
spara.printsteps = 50
spara.randomseed = 37
spara.restraintsfile = restraints
spara.savepdbsteps = 100
spara.shakeinitcrd = 1
spara.storetopn = 50
spara.temperaturegrps = ['ALL']
spara.temperatures = [1 * kbt]  # 单位kbt*T
spara.timestep = 0.002
spara.gamma = 0.5
spara.temperaturegrps = ['SELECTIONS', 'chain0', '0-50;', 'chain0', '51-100;', 'chain0', '101-229']
spara.temperatures = [2.4942, 4.9884, 7.4826]  # kbt*T

SCUBASD(spara, mpara, ipara, 500000)

###################Changeresidues###########################
pdbfile = 'sd1top_0.pdb'
outputfile = 'changesequence1.pdb'
cr = ChangeResidues(pdbfile)                        #读入pdb初始化
changedchainid = 0
#newsequence = 'TTCCPSIVHHHHHHVCRL'                 #Changesequence功能：序列长度要和pdb匹配
changesites = {0: 'H', 10: 'H', 15: 'C', 20: 'V'}
# cr.Changesequence(changedchainid,newsequence)     #替换序列：序列长度要和pdb匹配
cr.MakeLVG(changedchainid)                          #makeLVG
cr.Replaceresidues(changedchainid, changesites)     #替换残基
cr.writepdb(outputfile)                             #输出文件

#################loopsampler##########################
lspdb='start3.pdb'
parfile = 'ls.par'                              # 设定参数文件，如果没有，会自动生成模板参数文件并命名为此
ipara, sdpara, lpspara = loopsamplerinit(parfile, jobname='auto')
lpspara.pdbstart = lspdb                        # 输入的pdb文件 。必须
lpspara.reconstructnum = 2000                   # 进行蒙特卡洛loop采样的轮数，每一轮会重建1-3个loop区，进行一个短的恒温SD。重要
lpspara.findlooptimes = 30
lpspara.loopsearchmode = 'auto'                 # loop位置和长度的搜索模式，重要
lpspara.readloopnum = 2                         # 读取非冗余结构的数量
lpspara.readrmsd = '0.9'                        # 非冗余的RMSD标准
lpspara.outloopfile = 'outloops.dat'
loops = "A 73 81 0 A 104 112 0 A 128 138 0 "     # 需要采样的loop 。必须
setloop(lpspara, loops)

SCUBALoopSampler(parfile,lpspara, ipara, sdpara,loops)    # loop采样,并修改lsparfile

#autoreadloop(parfile)  #循环采样结束后读出采样的loop,只读默认值，可以使用SCUBALoopReader单独读取
###################loopreader##############
loopreadparfile = 'ls.par'
lpspara, sl = loopreaderinit(loopreadparfile)   #sl: saveloop
lpspara.readloopnum = 2                         #参数同上
lpspara.readrmsd = '0.5'

SCUBALoopReader(sl, lpspara)                    # 读loop构象替换到初始结构上，产生top_n.pdb和相关文件


#################designseq###################
pdbname='top_0.pdb'
outname ='SCUBAdesign'          #输出文件名：SCUBAdesign-00n.pdb
logfile='design.log'            #log文件
n= 2                            #设计序列产生数
resfile=''                      #resfile
parafile =''
aapropfile=''
#Designseq(pdbname, outname,logfile, n=1,vdw=1,div=1,nat=-100, resfile='', parafile='',aapropfile='')
Designseq(pdbname,outname,logfile,n)


#This program shows how to fix some residues in SD
from pySCUBA import chooseSDstrategy,runSDstrategy,Strategy,MolSystmPar,SDRunPar,InteractionPar


#############SD##############
sdpdb='build0.pdb'
# SSinfo='build0_final_ss.txt'
# restraints='restraints1.txt'
# writerestraints(sdpdb,SSinfo,restraints)  #根据pdb和二级结构信息（从sketch中生成）自动生成全部二级结构的rmsd，helix，scale，rg restraints，不包含contact restraints


usr_strategy,usr_para=chooseSDstrategy('default_strategy')      #选择SD策略，返还用户参数usr_para和用户策略usr_strategy
usr_para.MolSystmPar.PDBStart='build0.pdb'                      #选择输入的pdb文件
usr_para.SDRunPar.RestraintsFile='restraints1.txt'              #选择约束参数文件
usr_para.SDRunPar.AnnealingScheme='1.4 0.9 10000 2000 2000'     #模拟退火温度选项
usr_para.MolSystmPar.FixedResidues='chain0 0-50，60，59-73'     #选择固定的残基

runSDstrategy(usr_strategy, 50000)                              #运行选择的策略

# sdpar='sd_pdbbuild.par'                                       #或者直接给定参数文件
# SCUBASD(sdpar, 50000)                                         #按照参数文件运行SD

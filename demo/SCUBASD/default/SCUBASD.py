# This program show Global Stochastic Dynamic Simulated Annealing Sampling of Protein Structure
from pySCUBA import chooseSDstrategy,runSDstrategy,Strategy,MolSystmPar,SDRunPar,InteractionPar

kbt=2.4942
A2NM=0.1
#############SD##############
#输入文件
sdpdb='build0.pdb'                          #pdb必须
SSinfo='build0_final_ss.txt'                #二级结构文件，可选
restraints='restraints1.txt'                #约束文件
# writerestraints(sdpdb,SSinfo,restraints)  #自动生成全部二级结构的rmsd，helix，scale，rg restraints，不包含contact restraints


usr_strategy,usr_para=chooseSDstrategy('default_strategy')      #选择SD策略：default/build/relax/quench。这里选择的默认策略
usr_para.MolSystmPar.PDBStart=sdpdb                             #选择pdb
usr_para.SDRunPar.AnnealingScheme='2.4 0.9 10000 2000 2000'     #设置模拟退火温度选项，详见document
usr_para.SDRunPar.RestraintsFile=restraints                     #约束文件

runSDstrategy(usr_strategy, 500000)                             #runSD

# sdpar='sd_pdbbuild.par'
# SCUBASD(sdpar, 50000)                                         #也可以提供SDpar文件，按照所给的参数运行

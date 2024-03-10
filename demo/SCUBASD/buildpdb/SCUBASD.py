# This program show optimizing the structure generated from sketch
from pySCUBA import chooseSDstrategy,runSDstrategy,Strategy,MolSystmPar,SDRunPar,InteractionPar



#############SD##############
#输入文件
sdpdb='build0.pdb'                          #pdb必须
SSinfo='build0_final_ss.txt'                #二级结构文件，可选
restraints='restraints1.txt'                #约束文件
# writerestraints(sdpdb,SSinfo,restraints)  #自动生成全部二级结构的rmsd，helix，scale，rg restraints，不包含contact restraints


usr_strategy,usr_para=chooseSDstrategy('build_strategy')    #选择SD策略：default/build/relax/quench。这里选择的build策略，用来优化由Sketch生成的结构
usr_para.MolSystmPar.PDBStart=sdpdb                         #选择pdb
usr_para.SDRunPar.GroupTemperatures='0.4'                   #组分初始温度
usr_para.SDRunPar.RestraintsFile=restraints                 #约束文件

runSDstrategy(usr_strategy, 20000)                          #runSD

# sdpar='sd_pdbbuild.par'
# SCUBASD(sdpar, 50000)                                     #也可以提供SDpar文件，按照所给的参数运行

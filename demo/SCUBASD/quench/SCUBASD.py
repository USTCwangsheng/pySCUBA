# This program show quench strategy optimize structures down to very low temperatures,please set the number of total steps to 20000(the steps in one SD round you set)
from pySCUBA import chooseSDstrategy,runSDstrategy,Strategy,MolSystmPar,SDRunPar,InteractionPar

kbt=2.4942
A2NM=0.1
#############SD##############
#输入文件
sdpdb='build0.pdb'                          #pdb必须


usr_strategy,usr_para=chooseSDstrategy('quench_strategy')       #选择SD策略：default_strategy/build_strategy/relax_strategy/quench_strategy。这里选择的quench策略，会在一轮SD模拟中逐渐降低温度
usr_para.MolSystmPar.PDBStart=sdpdb
usr_para.SDRunPar.GroupTemperatures='0.4'
usr_para.SDRunPar.AnnealingScheme='0.5    0.1    20000 1000 5000'     #设置模拟退火温度选项，一轮20000的SD：温度从0.5kbt降到0.1kbt，其中1000步0.5kbt，5000步下降到0.1kbt，剩下14000步持续在0.1kbt


runSDstrategy(usr_strategy, 20000)                                  #runSD，选择quench_strategy时，步数应与一轮SD的步数相同。

# sdpar='sd_pdb.par'
# SCUBASD(sdpar, 50000)

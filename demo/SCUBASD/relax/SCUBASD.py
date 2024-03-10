# This program show constant low temperature SD simulation
from pySCUBA import chooseSDstrategy,runSDstrategy,Strategy,MolSystmPar,SDRunPar,InteractionPar

kbt=2.4942
A2NM=0.1
#############SD##############
sdpdb='build0.pdb'


usr_strategy,usr_para=chooseSDstrategy('relax_strategy')        #选择SD策略：default_strategy/build_strategy/relax_strategy/quench_strategy。这里选择的relax策略，持续低温模拟
usr_para.MolSystmPar.PDBStart=sdpdb
usr_para.SDRunPar.GroupTemperatures='0.4'                       #组分初始温度，relax_strategy没有模拟退火，温度保持在此基本不变


runSDstrategy(usr_strategy, 10000)                              #runSD

# sdpar='sd_pdbbuild.par'
# SCUBASD(sdpar, 50000)

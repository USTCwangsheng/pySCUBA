# This program show use sketch file generate pdb
from pySCUBA import Sketchbuildpdb

#################sketch生成初始结构#####################
sketchpar='sketch.par'#Parameter file, which is automatically generated if it does not exist
skecthfile='3pvh.txt'       #sketch coordinate file. Must exist.
linkoption=1                #Connect loop option. 0 is not connected and 1 is connected.
randomseed=15               
outputname='build'   #Output the file name of the pdb with this prefix.
gennumber=2                 #Number of structures generated
Sketchbuildpdb(sketchpar,skecthfile,linkoption,randomseed,outputname,gennumber)


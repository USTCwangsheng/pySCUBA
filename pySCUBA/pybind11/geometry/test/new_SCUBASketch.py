from lib_SCUBASketch import *
from pb_iblock import USRIntrctMol

####read parameters###

sketch=Sketch()

#option1:use para_file
sketch.read_parfile("timb.txt")
#option2:creat sketch
##sketch.add_element("E",[7,7],"N",[0,10.5,-4],[0,0,1])
sketch.add_element("E",[7,7],"N",[0,10.5,-4],[0,0,1])
sketch.add_element("H",[10,10],"C",[0,19,0],[0,0,-1])

####Reseed#######

RandomEngine.reseed(1)

####Loop#######
for i in range(2):
    imol=USRIntrctMol()
    modeler=MolModeler()
    modeler.settarget(imol)
    modeler.buildssele(*sketch.build)
    imol.usrwritepdb("ouput_{}.pdb".format(i))



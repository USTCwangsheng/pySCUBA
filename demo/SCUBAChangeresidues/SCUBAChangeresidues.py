from pySCUBA import ChangeResidues


#initialization
pdbfile='./3nir_pre.pdb'#pluto
cr=ChangeResidues(pdbfile)

#parameters
changedchainid=0
newsequence='TTCCPSIVHHHHHHVCRLPGHHHHHHATYTGCIIIPGATCPGDYAN'
changesites={0:'H',10:'H',15:'C',20:'V'}

#function
cr.Changesequence(changedchainid,newsequence)
cr.MakeLVG(changedchainid)
cr.Replaceresidues(changedchainid,changesites)

#output
outputfile='./changesequence.pdb'
cr.writepdb(outputfile)

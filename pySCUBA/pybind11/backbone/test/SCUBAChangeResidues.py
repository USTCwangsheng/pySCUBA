import sys
from pySCUBA import pb_iblock,pb_proteinrep

__all__=['ChangeResidues']
AAmap={
    "A":"ALA",
    "R":"ARG",
    "N":"ASN",
    "D":"ASP",
    "C":"CYS",
    "Q":"GLN",
    "E":"GLU",
    "G":"GLY",
    "H":"HIS",
    "I":"ILE",
    "L":"LEU",
    "K":"LYS",
    "M":"MET",
    "F":"PHE",
    "P":"PRO",
    "S":"SER",
    "T":"THR",
    "W":"TRP",
    "Y":"TYR",
    "V":"VAL",
}

def AAone2three(aa):
    aa=aa.upper()
    if aa not in AAmap.keys():
        return "XXX"
    else:
        return AAmap[aa]

def AAthree2one(aa):
    for one,three in AAmap.items():
        if aa.upper()==three:
            return one

    return "X"


class ChangeResidues:
    def __init__(self,startpdb):
        self.startpdb=startpdb

        #extract information about sequence and secondary structure
        rstates=pb_proteinrep.residuestates(self.startpdb)
        self.ssseq=[]
        self.seq=[]        
        for n in range(len(rstates)):
            seqidx=1
            ssseq=""
            seq=""
            for m in range(len(rstates[n])):
                seq+=pb_proteinrep.AminoAcidSeq.name2code(rstates[n][m].sidechainstate.resiudetype)
                ssseq+=rstates[n][m].backbonestate.SSState
                seqidx+=1
            self.seq.append(seq)
            self.ssseq.append(ssseq)

        self.modeler=pb_iblock.MolModeler()
        self.model=pb_proteinrep.AAConformersInModel()
        self.model.readpdbfile(self.startpdb)
        self.imol=pb_iblock.IntrctMol(self.model)
        self.modeler.settarget(self.imol)
        self.modeler.checkclash(True)

    def Changesequence(self,chainidx,newsequence):
        print("Changesequence")
        try:      
            for c in range(len(self.ssseq[chainidx])):
                AA=newsequence[c]
                self.modeler.changeresidue([chainidx,c],AAone2three(AA))
        except:
            print('FAILED: The {} contains error charactor or length'.format(chainidx))

    def MakeLVG(self,chainidx):
        print("MakeLVG")        
        for c in range(len(self.ssseq[chainidx])):
            if self.ssseq[chainidx][c]=="H": self.modeler.changeresidue([chainidx,c],"LEU")
            elif self.ssseq[chainidx][c]=="E": self.modeler.changeresidue([chainidx,c],"VAL")
            else: self.modeler.changeresidue([chainidx,c],"GLY")

    def Replaceresidues(self,chainidx,Changesites): #Changesites:dirtionary
        print("Replaceresidues")
        
        for n,aa in Changesites.items():
            if n>=len(self.ssseq[chainidx]):
                print("FAILED: The index {} of residue is ranged out!!!".format(n))
                return
            AA=aa[0]
            if AA not in AAmap.keys():
                print("FAILED: The amino acide type {} is not normal!!!".format(aa))
                return            
            self.modeler.changeresidue([chainidx,n],AAone2three(AA))

    def writepdb(self,outputfile):
        self.imol.writepdb(outputfile)


if __name__ == '__main__':
    #initialization
    pdbfile='/home/zhanglu/workspace/new/pybind11/backbone/test/3nir_pre.pdb'#pluto
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
    outputfile='/home/zhanglu/workspace/new/pybind11/backbone/test/changesequence.pdb' #pluto
    cr.writepdb(outputfile)




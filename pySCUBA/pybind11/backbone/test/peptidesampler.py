import argparse
from copy import deepcopy
from collections import namedtuple
import random
from math import sqrt
from pySCUBA import pb_backbone,pb_iblock,pb_sampling,pb_proteinrep,pb_designseq,pb_dstl,pb_geometry

__all__=['Protein','PeptideContainer','PepBuilder','PepSampler','PepSamplerPar','SDRunPar','InteractionPar']

class Protein:
    def __init__(self, pdbfile):
        self.pdbfile=pdbfile
        self.model=pb_proteinrep.AAConformersInModel()
        self.model.readpdbfile(self.pdbfile)
        self.imol=pb_iblock.IntrctMol(self.model)
        self.backbonesite=pb_backbone.readbackbonefrompdb(self.pdbfile)

    def nchains(self):
        return pb_iblock.nchains(self.imol)

    def nresidues(self,chainidx):
        return pb_iblock.nresidues(self.imol,chainidx)

    def getcacrd(self,chainidx,residueidx):
        return self.backbonesite[chainidx][residueidx].cacrd()

    def getcrds(self):
        natoms=self.imol.natoms()
        crds=pb_iblock.recollectcrds_all(self.imol)
        return crds

    def getcrds_double(self):
        self.imol.natoms()
        crds=pb_iblock.getcrds_all_d(self.imol)
        return crds

    def changecrds(self, crds):
        pb_iblock.changecrds(self.imol,crds)


    def changecrds_double(self, crds):
        assert(len(crds)%3==0)
        newcrds=[]
        for n in range(int(len(crds)/3)):
            newcrds.append(pb_geometry.XYZ(crds[n*3],crds[n*3+1],crds[n*3+2]))
        self.changecrds(newcrds)

    def addchain(self,blcks):
        pb_iblock.addchain(self.imol,blcks)
        self.imol.natoms()

    def writepdb(self,savedfile):
        self.imol.writepdb(savedfile)

    def specifyactiveblks(self,mainchainselectedblks,sidechainselectedblks):
        self.imol.specifyactiveblks(mainchainselectedblks,sidechainselectedblks)

    def setupintrctpara(self, ipara):
        pb_iblock.setintrctpara(self.imol, ipara)

    def getenergy(self):
        return pb_iblock.sum_energies(self.imol)


class PeptideContainer:
    def __init__(self,receptor,peplength,rmsdcutoff):
        self.receptor=receptor
        self.peplength=peplength
        self.rmsdcutoff=rmsdcutoff
        self.pepconfigs=[] #energies,totalenergy,contactrestEnergy,crds
    
    def savenonredundant(self,OneConfig):
        energysortid=-1
        ifsimilar=False
        for n in range(len(self.pepconfigs)):
            currentConfig=self.pepconfigs[n]
            if self.rmsd(currentConfig[-1],OneConfig[-1])<self.rmsdcutoff:
                if OneConfig[1]<currentConfig[1]:
                    self.pepconfigs[n]=OneConfig
                ifsimilar=True
            if OneConfig[1]<currentConfig[1]: energysortid=n
        
        if not ifsimilar:
            self.pepconfigs.insert(energysortid,OneConfig)

    def writepdb(self,savefile,oneconfig):
        xyzs=self.receptor.getcrds_double()
        for n in range(self.peplength):
            for i in range(4):
                for j in range(3):
                    xyzs[(-self.peplength+n)*12+i*3+j]=oneconfig[-1][n*12+i*3+j]
        
        self.receptor.changecrds_double(xyzs)
        self.receptor.writepdb(savefile)


    def rmsd(self,config1,config2): #config:  the list of xyz
        rmsdmin=10000.0
        minalignment=(self.peplength+1)/2
        residuedevs=[]
        for offset in range(int(self.peplength-minalignment+1)):
            for i in range(int(self.peplength-offset)):
                diff2=0.0
                for m in range(4):
                    for z in range(3):
                        # dx=config1[4*i+m]-config2[4*i+m+4*offset]
                        dx=config1[12*i+3*m+z]-config2[12*i+3*m+z+12*offset]
                        diff2+=dx**2
                residuedevs.append(diff2)

            residuedevs.sort()
            rmsd=0.0
            for i in range(int(minalignment)):
                rmsd+=residuedevs[i]

            rmsd=sqrt(rmsd/(4.0*minalignment))
            if rmsd<rmsdmin: rmsdmin=rmsd

        return rmsd

    def writepdb_topn(self,n,rmsdcutoff): #TODO
        currentpepconfigs=deepcopy(self.pepconfigs)
        delidxlist=[]
        for i in range(len(self.pepconfigs)-1):
            for m in range(i+1,len(self.pepconfigs)):
                if self.rmsd(self.pepconfigs[i][-1],self.pepconfigs[m][-1]) < rmsdcutoff:
                    if self.pepconfigs[i][1]<self.pepconfigs[m][1]:
                        delidxlist.append(m)
                        # del currentpepconfigs[m]

        counter=0
        if len(currentpepconfigs)<n: n=len(currentpepconfigs)
        for m in range(len(currentpepconfigs)):
            if m in delidxlist: continue
            counter+=1
            self.writepdb(str(m)+".pdb",currentpepconfigs[m])
            if counter==n: break


PepSamplerPar=namedtuple("PepSamplerPar",field_names=[
    "RandomSeed","MaxSDSteps","PrintParameters","ContactRestraint","ENEAverageDecay","ENEVarStop",
    "JobName","Verbose","SaveConfigFile","LigPepLength","InterfaceResidues","PDBStart","RMSDCutOff"
],defaults=[
    "1357","5000","0","5 500.0","0.98","1",
    "lig0","0","savedligconfig.dat",10,"","","0.0"
])

SDRunPar=namedtuple("SDRunPar",field_names=[
    "DOAnnealing","AnnealingScheme","DOShake","PrintParameters","PrintSteps","RandomSeed",
    "RecalcNeighborListSteps","RecalcSSSteps","SavePDBSteps","GAMMA","TimeStep","GroupTemperatures",
    "OutPDBFile","TemperatureGroups"
],defaults=[
    0,[0.5,0.1,5000,1000,3000],1,0,50,513,
    100,500,10000,0.5,0.002,[0.1],
    "sdtest_out.pdb",["ALL"]
])

InteractionPar=namedtuple("InteractionPar",field_names=[
    'CovalentWeight','MCStericWeight','PhiPsiWeight','LocalStrWeight','RotamerWeight',
    'SitePairWeight','SCPackingWeight','LocalHBWeight'
],defaults=[
    1.0,1.4,2.0,0.5,2.4,0.32,3.1,0.6
])

def peptideparas():
    return PepSamplerPar(),SDRunPar(),InteractionPar()

class SDRunPara:
    def __init__(self,SDRunPara) -> pb_sampling.SDRunPara():
        self.sdpara=pb_sampling.SDRunPara()
        self.sdpara.doannealing=SDRunPara.DOAnnealing
        self.sdpara.annealingscheme=SDRunPara.AnnealingScheme
        self.sdpara.doshake=SDRunPara.DOShake
        self.sdpara.printsteps=SDRunPara.PrintSteps
        self.sdpara.randomseed=SDRunPara.RandomSeed
        self.sdpara.savepdbsteps=SDRunPara.SavePDBSteps
        self.sdpara.gamma=SDRunPara.GAMMA
        self.sdpara.timestep=SDRunPara.TimeStep
        self.sdpara.temperaturegrps=SDRunPara.TemperatureGroups
        self.sdpara.temperatures=SDRunPara.GroupTemperatures
        self.sdpara.outpdbfile=SDRunPara.OutPDBFile
        
        
class IntrctPara:
    def __init__(self,InteractionPar) -> pb_iblock.IntrctPara():
        self.ipara=pb_iblock.IntrctPara()
        self.ipara.weight_coval=InteractionPar.CovalentWeight
        self.ipara.weight_steric=InteractionPar.MCStericWeight
        self.ipara.weight_phipsi=InteractionPar.PhiPsiWeight
        self.ipara.weight_localstr=InteractionPar.LocalStrWeight
        self.ipara.weight_rotamer=InteractionPar.RotamerWeight
        self.ipara.weight_sitepair=InteractionPar.SitePairWeight
        self.ipara.weight_scpacking=InteractionPar.SCPackingWeight
        self.ipara.weight_localhb=InteractionPar.LocalHBWeight
        

class PepBuilder:
    def __init__(self,receptor,peptidelen,ifresiduesidx):
        self.receptor=receptor
        self.peptidelen=peptidelen
        self.backbonebuilder=pb_backbone.BackBoneBuilder()
        self.ifresiduesidx=[]
        ifreslist=ifresiduesidx.split()
        for n in range(1,len(ifreslist)):
            self.ifresiduesidx.append([int(ifreslist[0][-1]),n])
        totalres=0
        self.receptorcenter=pb_geometry.XYZ(0.0,0.0,0.0)
        for nchain in range(self.receptor.nchains()):
            for nres in range(self.receptor.nresidues(nchain)):
                self.receptorcenter+=self.receptor.getcacrd(nchain,nres)
                totalres+=1
        self.receptorcenter/=totalres

        self.interfacerescacrds=[]
        for idx1 in range(len(self.ifresiduesidx)-1):
            for idx2 in range(idx1+1,len(self.ifresiduesidx)):
                firstcacrd=self.receptor.getcacrd(self.ifresiduesidx[idx1][0],self.ifresiduesidx[idx1][1])
                secondcacrd=self.receptor.getcacrd(self.ifresiduesidx[idx2][0],self.ifresiduesidx[idx2][1])
                dis=firstcacrd.distance(secondcacrd)
                self.interfacerescacrds.append([self.ifresiduesidx[idx1][1],self.ifresiduesidx[idx2][1],firstcacrd,secondcacrd,dis])

        assert(len(self.interfacerescacrds)>=2)

    def genpepcrds(self):
        center,direction=self.randompepcenter()
        # print('center crd:'+center.__repr__())
        # print('direction crd:'+direction.__repr__())
        bs0=pb_backbone.BackBoneSite()
        pb_backbone.genbackbonesite(False,-120.0,120.0,bs0) #judge bs1
        helixregion=[]
        strandregion=[]
        cissites=[]          
        bss=pb_backbone.buildforwardbackbone(self.peptidelen,bs0,helixregion,strandregion,cissites)
        xyz=pb_geometry.XYZ(0.0,0.0,0.0)
        newbss=pb_backbone.movechainto(xyz,direction,True,bss)
        xyzs=[]
        for bs in newbss:
            for i in range(4):
                xyzs.append(bs.getcrd(int(pb_backbone.BackBoneSite.Field.NCRD)+3*i))
        oldc=pb_geometry.XYZ(0.0,0.0,0.0)
        for r in xyzs: oldc=oldc+r
        oldc=center-(1/(4.0*self.peptidelen))*oldc
        for r in xyzs: r=0.1*(r+oldc)
        return xyzs

    def genpepcrd(self,receptor):
        xyzs=self.genpepcrds()
        crds=receptor.getcrds()
        for n in range(-(len(xyzs)),0):
            crds[n]=xyzs[n+len(xyzs)]/10.0
        receptor.changecrds(crds)

    def addpeptoreceptor(self,receptor,gencrd):
        modeler=pb_iblock.MolModeler()
        modeler.settarget(receptor.imol)
        blcks=[]
        bs=pb_backbone.BackBoneSite()
        for i in range(self.peptidelen):
            blck=pb_iblock.IntrctBlck()
            blck=modeler.mkfloatingblck([-10,-10],'GLY',bs)
            blcks.append(blck)

        receptor.addchain(blcks)
        
        if gencrd:
            # cid=receptor.nchains()-1
            self.genpepcrd(receptor)

    def randomxyz(self,ra):
        get=False
        while not get:        
            norm=0.0
            p=[]
            for d in range(3):
                x=random.uniform(0.0,1.0)  #uniform random number between 0 to 1.0
                x=x*2.0-1.0
                p.append(x)
                norm +=x*x                
            if norm >1.0: continue
            get=True
            x=p[0]*ra
            y=p[1]*ra
            z=p[2]*ra
        point=pb_geometry.XYZ(x,y,z)
        return point

    def randompepcenter(self):
        idx=random.randint(0,len(self.interfacerescacrds)-1)
        judgenum=random.randrange(0.0,1.0)
        rm=(self.interfacerescacrds[idx][2]+self.interfacerescacrds[idx][3])/2
        rd=rm-self.receptorcenter
        rd=(1.0/sqrt(rd.squarednorm()))*rd
        center=rm+rd+self.randomxyz(0.1)
        if(judgenum<0.5): direction=self.interfacerescacrds[idx][2]-self.interfacerescacrds[idx][3]
        else: direction=self.interfacerescacrds[idx][3]-self.interfacerescacrds[idx][2]
        direction=direction+self.randomxyz(0.1*(sqrt(direction.squarednorm())))
        
        return center, direction

class PepSampler:
    def __init__(self, receptor, pepsamplerpar, sdpar, ipara):
        # self.receptor=receptor #Protein object        
        self.pepsamplerpar=pepsamplerpar
        # self.sdpar=sdpar
        self.pepbuilder=PepBuilder(receptor,self.pepsamplerpar.LigPepLength,pepsamplerpar.InterfaceResidues)
        self.pepbuilder.addpeptoreceptor(receptor,True)
        self.receptor=receptor        
        self.optimizer=pb_sampling.MolSDRun()
        self.newsdpar=SDRunPara(sdpar)
        self.ipara=IntrctPara(ipara)
        self.receptor.setupintrctpara(self.ipara.ipara)
        self.setupsdrun(self.newsdpar.sdpara)
        self.pepcontainers=PeptideContainer(self.receptor,self.pepsamplerpar.LigPepLength,float(self.pepsamplerpar.RMSDCutOff))


    def buildandoptimize(self):
        self.rebuildpeptide()
        initcrd=self.receptor.getcrds_double()
        self.optimizer.reinitstate(initcrd,0)
        callback=pb_sampling.SDCallBack_default(self.newsdpar.sdpara)
        self.optimizer.runsteps(int(self.pepsamplerpar.MaxSDSteps),callback)


    def saveconfig(self,savefile):
        etotlist=self.receptor.getenergy()
        etot=etotlist[0]
        with open(savefile,'a') as wdata:
            wdata.write("{:.2f}\n".format(float(etot)))
            for e in etotlist[1:]: wdata.write("{:.2f} ".format(float(e)))
            wdata.write('\n')
            for eneres in pb_sampling.restraints_energies(self.optimizer): wdata.write("{:.2f} ".format(float(eneres)))
            wdata.write('\n')
            xyzs=self.receptor.getcrds_double()
            pepcrds=[]
             
            for nbegin in range(self.pepsamplerpar.LigPepLength,0,-1):
                for i in range(4):
                    wdata.write("{:.2f} {:.2f} {:.2f}\n".format(float(xyzs[len(xyzs)-nbegin*12+i*3]),float(xyzs[len(xyzs)-nbegin*12+i*3+1]),float(xyzs[len(xyzs)-nbegin*12+i*3+2])))
                    pepcrds.append(xyzs[len(xyzs)-nbegin*12+i*3])
                    pepcrds.append(xyzs[len(xyzs)-nbegin*12+i*3+1])
                    pepcrds.append(xyzs[len(xyzs)-nbegin*12+i*3+2])
            
            #pepcontainers
            self.pepcontainers.savenonredundant([etotlist[1:],etot,pb_sampling.restraints_energies(self.optimizer)[0],pepcrds]) #energies,totalenergy,contactrestEnergy,crds

    def writeContainers(self,n,rmsdcutoff):
        self.pepcontainers.writepdb_topn(n,rmsdcutoff)

    def setupsdrun(self,sdpara):
        pepcid=self.receptor.nchains()-1
        ifreslist=self.pepsamplerpar.InterfaceResidues.split()
        resblcks=[]
        resblcks_lig=[]
        for n in range(1,len(ifreslist)):
            resblcks.append(int(ifreslist[n])) #[int(ifreslist[0][-1]),resblcks]
        for n in range(self.receptor.nresidues(pepcid)):
            resblcks_lig.append(n)

        iresiduesreceptor=[[int(ifreslist[0][-1]),set(resblcks)]]
        iresiduesligand=[[int(pepcid),set(resblcks_lig)]]

        self.receptor.specifyactiveblks(iresiduesligand,[])
        self.optimizer.init(sdpara,self.receptor.imol,False)

        #contact restraint
        contactres=pb_sampling.ContactRestraint()        
        nc0_i=5.0
        kres_i= 500.0
        resgdmin_i= 0.5
        resgdsmall_i= 0.8
        resgdoff_i= 2.0
        contactres=pb_sampling.initContactRestraint(self.receptor.imol, iresiduesreceptor, iresiduesligand, nc0_i, kres_i, resgdmin_i, resgdsmall_i, resgdoff_i)
        self.optimizer.addrestraint_contact(contactres)


    def rebuildpeptide(self):
        self.pepbuilder.genpepcrd(self.receptor)


def args():
    parser=argparse.ArgumentParser(argument_default=argparse.SUPPRESS)
    parser.add_argument('--nconfig',type=int,help='the number of sampling proteins')
    parser.add_argument('--rmsdcutoff',type=float,help='the rmsd cutoff among those stored peptides')
    parser.add_argument('--outputnum',type=int,help='the output number of peptide')

    args=vars(parser.parse_args())
    return args

def setseed(seed):
    random.seed(seed)

if __name__ == '__main__':
    para=args()    

    #protein object
    pdbfile='/home/zhanglu/workspace/new/pybind11/backbone/test/1awr_preprocessed.pdb' #pluto
    protein=Protein(pdbfile)

    interfaceresidues='chain0 71 72 101 102 62 100'
    pepsamplerpar=PepSamplerPar(PDBStart=pdbfile,InterfaceResidues=interfaceresidues)

    sdpar=SDRunPar()
    ipar=InteractionPar()
    setseed(int(sdpar.RandomSeed))

    #test pepbuilder program
    # savefile='/home/zhanglu/workspace/new/pybind11/backbone/test/addpeptide.pdb'
    # pepbuilder=PepBuilder(protein,pepsamplerpar.LigPepLength,pepsamplerpar.InterfaceResidues)
    # pepbuilder.addpeptoreceptor(protein,true)
    # protein.writepdb(savefile)


    #test pepsampler program
    savefile='/home/zhanglu/workspace/new/pybind11/backbone/test/pepsavedconfigs'
    pepsampler=PepSampler(protein,pepsamplerpar,sdpar,ipar)
    for n in range(para['nconfig']):
        pepsampler.buildandoptimize()
        pepsampler.saveconfig(savefile)
    pepsampler.writeContainers(para['outputnum'],para['rmsdcutoff'])
    print("Done!")
    
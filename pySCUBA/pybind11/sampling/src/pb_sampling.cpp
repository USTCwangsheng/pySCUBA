/*
 * @Author: zjw
 * @2022-03-31
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>
#include <stdlib.h>
#include <sampling/groupdistance.h>
#include <sampling/ligandpep.h>
#include <sampling/ligpepsampler.h>
#include <sampling/loopsampler.h>
#include <sampling/molsdrun.h>
#include <sampling/sdrunpara.h>
#include <sampling/runsddefaults.h>
#include <sampling/restraints.h>
#include <sampling/stochasticdynamics.h>
#include <sampling/pocketdocker.h>
#include <sampling/screen.h>
#include <sampling/savedloops.h>
#include <sampling/test/SCUBALoopSampler.cpp>
//#include <sampling/test/SCUBALoopReader.cpp>
#include <geometry/rectgrid.h>
#include <geometry/sphere.h>
//#include <sampling/SCUBALoopSampler.cpp>
namespace py = pybind11;
namespace spl = NSPsampling;
using namespace NSPintrct;
typedef std::pair<int, NSPgeometry::XYZ> DvDxi;
typedef NSPintrct::MolModeler::SiteIdx SiteIdx;
typedef std::pair<char,std::string> PdbPosiIdx;
typedef NSPintrct::MolModeler::LoopConfMaker ConfMaker;
typedef std::pair<SiteIdx,SiteIdx>FlankingSite;

//SCUBALoopSampler ws
void pb_backupdata(NSPsampling::LoopSampler& sampler) {
    sampler.imol()->getmolsystm().backupdata();
}
void pb_dropbackdata(NSPsampling::LoopSampler& sampler) {
    sampler.imol()->getmolsystm().dropbackdata();
}

void pb_readloop(NSPsampling::LpSamplerPara& lpspara, std::vector<std::string>& loops) {
    std::vector<std::pair<PdbPosiIdx, PdbPosiIdx>> emp_flankingpositions;
    std::vector<int> emp_looplengths;
    //LpSamplerPara lpspara;
    assert(loops.size() % 4 == 0);
    int idx = 0;
    lpspara.nloops = loops.size() / 4;
    lpspara.flankingpositions= emp_flankingpositions;
    lpspara.looplengths= emp_looplengths;
    for (int i = 0; i < lpspara.nloops; i++) {
        char cid = loops[idx++][0];
        if (cid == '0') cid = ' ';
        std::string left = loops[idx++];
        std::string right = loops[idx++];
        int length = std::stoi(loops[idx++]);
        lpspara.flankingpositions.push_back({ {cid,left},{cid,right} });
        std::cout << "Loop " << i << " chain " << cid << " " << left << "- " << right << " length " << length << std::endl;
        lpspara.looplengths.push_back(length);
    }
}

LpSamplerPara makeloopreaderpar(SavedLoops& sl, std::string& parfile) {
    ControlFile cf;
    cf.readfile(parfile);
    LpSamplerPara lpspara = makelpsamplerparam(std::string(), cf);
    sl.init(lpspara);
    return lpspara;
}
void writestartpdb(SavedLoops& sl, std::string& startpdbconfig) {
    std::ofstream ofs(startpdbconfig);
    sl.writepdb(ofs, sl.startconfig, true);
}


void readloopout(SavedLoops& sl, LpSamplerPara& lpspara) {
    int ntop = 100;
    if (lpspara.readloopnum > 0)
        ntop = lpspara.readloopnum;
    if (lpspara.readrmsd == "auto")
    {
        std::cout << "Will use " << sl.rmsdcut() * 10 << " Angstrom for RMSD_Cutoff" << std::endl;
        sl.savenonredundant("nonreduntloops.dat", "nr_rmsdene.dat", sl.rmsdcut());
        sl.writepdb_topn(ntop, sl.rmsdcut());
    }
    else if (lpspara.readrmsd == "min")
    {
        double nla = (double)sl.nloopatoms;
        double denominator = nla / sl.min_nloopatom;
        std::cout << "Will use " << sl.rmsdcut_de(denominator) * 10 << " Angstrom for RMSD_Cutoff" << std::endl;
        sl.savenonredundant("nonreduntloops.dat", "nr_rmsdene.dat", sl.rmsdcut_de(denominator));
        sl.writepdb_topn(ntop, sl.rmsdcut_de(denominator));
    }
    else
    {
        double rmsdcut = 0.0;
        //std::cout << "lpspara.readrmsd " << lpspara.readrmsd << std::endl;
        //assert(rmsdcut = std::stod(lpspara.readrmsd));
        rmsdcut = atof(lpspara.readrmsd.c_str());
        std::cout << "Will use " << rmsdcut << " Angstrom for RMSD_Cutoff" << std::endl;
        sl.savenonredundant("nonreduntloops.dat", "nr_rmsdene.dat", rmsdcut / 10);
        sl.writepdb_topn(ntop, sl.rmsdcut_de(rmsdcut / 10));
    }
}


struct Returnoptloop {
    std::vector<NSPgeometry::XYZ> optloops_;
    double ene_;
    std::vector<double> energies_;
    //std::vector<std::vector<double>> loopenergies_;
    std::vector<std::array<double, NSPintrct::IntrctBlck::ENESIZE>> loopenergies_;
};

Returnoptloop optimizeloops(spl::LoopSampler lsr, double& ene, std::array<double, NSPintrct::IntrctBlck::ENESIZE>& energies) {
    Returnoptloop optl;
    optl.optloops_ = *lsr.optimizeloops(ene, energies);
    optl.ene_ = ene;
    std::vector<double> vx(energies.begin(), energies.end());
    optl.energies_ = vx;
    return optl;
    //return *lsr.optimizeloops(ene,energies);
}
std::vector<NSPgeometry::XYZ> loopatomcrds(spl::LoopSampler lsr) {
    return *lsr.loopatomcrds();
}
std::vector<NSPgeometry::XYZ> regenloops(spl::LoopSampler lsr,double &ene,std::array<double,NSPintrct::IntrctBlck::ENESIZE> &energies){
    return *lsr.regenloops(ene,energies);
}
Returnoptloop regenuseloops(spl::LoopSampler lsr,double &ene,std::array<double,NSPintrct::IntrctBlck::ENESIZE> &energies){
    Returnoptloop optl;
    optl.optloops_ = *lsr.regenuseloops(ene, energies);
    optl.ene_ = ene;
    std::vector<double> vx(energies.begin(), energies.end());
    optl.energies_ = vx;
    return optl;
    //return *lsr.regenuseloops(ene,energies);
}
Returnoptloop loopenes(spl::LoopSampler lsr, std::vector<std::array<double, NSPintrct::IntrctBlck::ENESIZE>>& loopenergies) {
    Returnoptloop optl;
    optl.energies_ = lsr.loopenes(loopenergies);

    optl.loopenergies_ = loopenergies;
    return optl;
    //return *lsr.regenuseloops(ene,energies);
}
NSPintrct::IntrctMol imol(spl::LoopSampler lsr) {
    return *lsr.imol();
}
//molsdrun.h
spl::StochasticDynamics sd(spl::MolSDRun msr) {
		return *msr.sd();
}
NSPintrct::IntrctMol imol(spl::MolSDRun msr) {
		return *msr.imol();
	}
spl::ShakeBonds shakebds(spl::MolSDRun msr) {
		return *msr.shakebds();
}

void scubasd(spl::SDRunPara spara, NSPintrct::MolSystmPara mpara, NSPintrct::IntrctPara ipara, int nsteps) {
    spl::SDCallBack_default callback(spara);
    spl::MolSDRun sdrun = spl::makemolsdrun(spara, mpara, ipara);
    sdrun.runsteps(nsteps, callback);
}

void SCUBASDparmod(std::string& parfile, int& nsteps) {
    std::string jobname("auto");
    NSPintrct::MolSystmPara mpara;
    NSPintrct::IntrctPara ipara;
    NSPsampling::SDRunPara spara;
    readparameters(parfile, mpara, ipara, spara, jobname);
    NSPsampling::SDCallBack_default callback(spara);
    NSPsampling::MolSDRun sdrun = makemolsdrun(spara, mpara, ipara);
    sdrun.runsteps(nsteps, callback);
}
void realoopfrompar(std::string& parfile) {
    SavedLoops sl;
    ControlFile cf;
    cf.readfile(parfile);
    LpSamplerPara lpspara = makelpsamplerparam(std::string(), cf);
    sl.init(lpspara);
    sl.writermsdene("rmsdene.dat");
    std::string startpdbconfig = "nativequenched.pdb";
    writestartpdb(sl, startpdbconfig);
    readloopout(sl, lpspara);
}

void loopsamplefrompar(std::string& parfile) {
    std::string jobname = "auto";
    IntrctPara ipara;
    SDRunPara sdpara;
    LpSamplerPara lpspara;
    readpara_lpsampler(parfile, lpspara, ipara, sdpara, jobname);
    searchloopandsample(lpspara, ipara, sdpara);
    autoreadloop(parfile);
}



PYBIND11_MAKE_OPAQUE(std::vector<NSPgeometry::XYZ>);
template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

PYBIND11_MODULE(pb_sampling, m)
{   
    //additional function
    m.def("scubasd", &scubasd);
    m.def("SCUBASDparmod", &SCUBASDparmod);
    m.def("realoopfrompar", &realoopfrompar);
    m.def("loopsamplefrompar", &loopsamplefrompar);
    m.def("printenergies", &printenergies);
    m.def("writeloops", &writeloops);
    m.def("backupdata", &pb_backupdata);
    m.def("dropbackdata", &pb_dropbackdata);
    m.def("writeloops", &writeloops);
    m.def("baklspar", &baklspar);
    m.def("autoreadloop", &autoreadloop);
    m.def("optimizeloops", &optimizeloops);
    m.def("loopenes", &loopenes);
    m.def("regenuseloops", &regenuseloops);
    m.def("readloop", &pb_readloop);
    m.def("makeloopreaderpar", &makeloopreaderpar);
    m.def("writestartpdb", &writestartpdb);
    m.def("readloopout", &readloopout);
    m.def("searchloopandsample", &searchloopandsample);
    //m.def("ofstream", &pb_ofstream);
    //m.def("writepdb", &pb_writepdb);

    py::class_<std::ofstream>(m, "ofstream")
        .def(py::init<std::string&>());

    py::class_<Returnoptloop>(m, "Returnoptloop")
        .def(py::init<>())
        .def_readwrite("optloops_", &Returnoptloop::optloops_)
        .def_readwrite("ene_", &Returnoptloop::ene_)
        .def_readwrite("loopenergies_", &Returnoptloop::loopenergies_)
        .def_readwrite("energies_", &Returnoptloop::energies_);

    //groupdistance.h
	py::class_<spl::GroupDistance>(m, "GroupDistance")
        .def(py::init<>())
        .def_readwrite("grps", &spl::GroupDistance::grps)
		.def("distance", &spl::GroupDistance::distance);
    
    py::class_<spl::GroupContact>(m, "GroupContact")
        .def(py::init<>())
        .def_readwrite("grpd", &spl::GroupContact::grpd)
        .def_readwrite("gdmin", &spl::GroupContact::gdmin)
        .def_readwrite("gdsmall", &spl::GroupContact::gdsmall)
        .def_readwrite("gdoff", &spl::GroupContact::gdoff)
        .def_readwrite("thetasmall", &spl::GroupContact::thetasmall)
		.def("theta", overload_cast_<const std::vector<NSPgeometry::XYZ> &,std::vector<DvDxi> *>()(&spl::GroupContact::theta,py::const_))
        .def("theta", overload_cast_<double,double *>()(&spl::GroupContact::theta,py::const_));
    
    //restraints.h zl
    py::class_<spl::ContactRestraint>(m, "ContactRestraint")
        .def(py::init<>())
    ;

    m.def("initContactRestraint",[](NSPintrct::IntrctMol & imol,
        std::vector<std::pair<int,std::set<int>>> iresiduesreceptor,
        std::vector<std::pair<int,std::set<int>>> iresiduesligand,
        double nc0_i,double kres_i, double resgdmin_i, double resgdsmall_i, double resgdoff_i){
            NSPintrct::BlkSelector iresidues_receptor,iresidues_ligand;
            for(int n=0;n<iresiduesreceptor.size();n++){
                iresidues_receptor.selectedblks[iresiduesreceptor[n].first]=iresiduesreceptor[n].second;
            }
            for(int n=0;n<iresiduesligand.size();n++){
                iresidues_ligand.selectedblks[iresiduesligand[n].first]=iresiduesligand[n].second;
            }
            return ContactRestraint(imol, iresidues_receptor, iresidues_ligand, nc0_i, kres_i, resgdmin_i, resgdsmall_i, resgdoff_i);
    });

    //ligandpep.h
    py::class_<spl::LigandPep>(m, "LigandPep")
        .def(py::init<>())
        //.def(py::init<const NSPintrct::IntrctMol &,const std::vector<NSPdstl::Idx2D> &>())
        .def("setup", &spl::LigandPep::setup)
        .def("genligandcrds", &spl::LigandPep::genligandcrds)
        .def("randomligandcenter", &spl::LigandPep::randomligandcenter)
        .def("addligandpep", &spl::LigandPep::addligandpep)
        .def("genligandpepcrd", &spl::LigandPep::genligandpepcrd);
    py::class_<spl::LigPConfigRMSD>(m, "LigPConfigRMSD")
        .def(py::init<int,int>())
        .def_readwrite("offsetmin", &spl::LigPConfigRMSD::offsetmin)
        .def("__call__", &spl::LigPConfigRMSD::operator());
    
    //ligpepsampler.h
    py::class_<spl::LigPSamplerPara>(m, "LigPSamplerPara")
        .def(py::init<const std::vector<std::string> &>())
        .def_readwrite("jobname", &spl::LigPSamplerPara::jobname)
        .def_readwrite("pdbstart", &spl::LigPSamplerPara::pdbstart)
        .def_readwrite("ifresidues", &spl::LigPSamplerPara::ifresidues)
        .def_readwrite("saveconfigfile", &spl::LigPSamplerPara::saveconfigfile)
        .def_readwrite("ligpeplength", &spl::LigPSamplerPara::ligpeplength)
        .def_readwrite("ncontacts0", &spl::LigPSamplerPara::ncontacts0)
        .def_readwrite("ligpeplength", &spl::LigPSamplerPara::ligpeplength)
        .def_readwrite("kcontactres", &spl::LigPSamplerPara::kcontactres)
        .def_readwrite("randomseed", &spl::LigPSamplerPara::randomseed)
        .def_readwrite("sdprintsteps", &spl::LigPSamplerPara::sdprintsteps)
        .def_readwrite("verbose", &spl::LigPSamplerPara::verbose)
        .def_readwrite("maxsdsteps", &spl::LigPSamplerPara::maxsdsteps)
        .def_readwrite("enedecay", &spl::LigPSamplerPara::enedecay)
        .def_readwrite("sdenevarcut", &spl::LigPSamplerPara::sdenevarcut)
        .def("mksdjudger", &spl::LigPSamplerPara::mksdjudger);
    m.def("makeligpsamplerparam",&spl::makeligpsamplerparam);
    m.def("readpara_ligpepsampler",&spl::readpara_ligpepsampler);
    py::class_<spl::LigPepSampler>(m, "LigPepSampler")
        .def(py::init<>())
        .def("setupmol",&spl::LigPepSampler::setupmol)
        .def("imol",&spl::LigPepSampler::imol)
        .def("setupsdrun",&spl::LigPepSampler::setupsdrun)
        .def("rebuildligand",&spl::LigPepSampler::rebuildligand)
        .def("buildandoptimize",&spl::LigPepSampler::buildandoptimize)
        .def("saveconfig",&spl::LigPepSampler::saveconfig)
        .def("mypara_mutable", overload_cast_<>()(&spl::LigPepSampler::mypara))
        .def("mypara_const", overload_cast_<>()(&spl::LigPepSampler::mypara, py::const_));
    m.def("mkligpepsampler",&spl::mkligpepsampler);

    //loopsampler.h
    //m.def("printenergies", &printenergies);
    m.def("writeloops", &writeloops);


    py::class_<spl::LpSamplerPara>(m, "LpSamplerPara")
        .def(py::init<>())
        .def(py::init<const std::vector<std::string> &>())
        .def_readwrite("flankingpositions", &spl::LpSamplerPara::flankingpositions)
        .def_readwrite("nloops", &spl::LpSamplerPara::nloops)
        .def_readwrite("flankingsites", &spl::LpSamplerPara::flankingsites)
        .def_readwrite("looplengths", &spl::LpSamplerPara::looplengths)
        .def_readwrite("flankingsites", &spl::LpSamplerPara::flankingsites)
        .def_readwrite("jobname", &spl::LpSamplerPara::jobname)
        .def_readwrite("flankingsites", &spl::LpSamplerPara::flankingsites)
        .def_readwrite("pdbstart", &spl::LpSamplerPara::pdbstart)
        .def_readwrite("outloopfile", &spl::LpSamplerPara::outloopfile)
        .def_readwrite("mcrandomseed", &spl::LpSamplerPara::mcrandomseed)
        .def_readwrite("sdprintsteps", &spl::LpSamplerPara::sdprintsteps)
        .def_readwrite("verbose", &spl::LpSamplerPara::verbose)
        .def_readwrite("mckbt", &spl::LpSamplerPara::mckbt)
        .def_readwrite("maxcycles", &spl::LpSamplerPara::maxcycles)
        .def_readwrite("maxsdsteps", &spl::LpSamplerPara::maxsdsteps)
        .def_readwrite("enedecay", &spl::LpSamplerPara::enedecay)
        .def_readwrite("sdenevarcut", &spl::LpSamplerPara::sdenevarcut)
        .def_readwrite("pmax", &spl::LpSamplerPara::pmax)
        .def_readwrite("unchangedcut", &spl::LpSamplerPara::unchangedcut)
        .def_readwrite("poolsizecut", &spl::LpSamplerPara::poolsizecut)
        .def_readwrite("reconstructnum", &spl::LpSamplerPara::reconstructnum)
        .def_readwrite("helixinloop", &spl::LpSamplerPara::helixinloop)
        .def_readwrite("readloopnum", &spl::LpSamplerPara::readloopnum)
        .def_readwrite("readrmsd", &spl::LpSamplerPara::readrmsd)
        .def_readwrite("phipsirmsd", &spl::LpSamplerPara::phipsirmsd)
        .def_readwrite("loopsearchmode", &spl::LpSamplerPara::loopsearchmode)
        .def_readwrite("findlooptimes", &spl::LpSamplerPara::findlooptimes)//ws
        .def("mkmcjudger",&spl::LpSamplerPara::mkmcjudger)
        .def("mksdjudger",&spl::LpSamplerPara::mksdjudger);
    py::class_<spl::LoopPool>(m, "LoopPool")
        .def_readwrite("loopconfs", &spl::LoopPool::loopconfs)
        .def_readwrite("p_frompool", &spl::LoopPool::p_frompool)
        .def_readwrite("unchanged", &spl::LoopPool::unchanged)
        .def_readwrite("P_max", &spl::LoopPool::P_max)
        .def_readwrite("Unchanged_th", &spl::LoopPool::Unchanged_th)
        .def_readwrite("Size_th", &spl::LoopPool::Size_th);
    py::bind_vector<std::vector<NSPgeometry::XYZ>>(m, "VectorXYZ");
    py::class_<spl::LoopSampler>(m, "LoopSampler")
        .def(py::init<>())
        .def("setupmol",&spl::LoopSampler::setupmol)
        .def("setupmolfind",&spl::LoopSampler::setupmolfind)
        .def("setupsdrun",&spl::LoopSampler::setupsdrun)
        .def("loopblks",&spl::LoopSampler::loopblks)
        .def("findloops",&spl::LoopSampler::findloops)
        .def("findshortloops", &spl::LoopSampler::findshortloops)//ws
        .def("fixshortloops", &spl::LoopSampler::fixshortloops)//
        .def("fixloops", &spl::LoopSampler::fixloops)//
        .def("saveloops",&spl::LoopSampler::saveloops)
        .def("rebuildloops",&spl::LoopSampler::rebuildloops)
        .def("optimizeloops", overload_cast_<>()(&spl::LoopSampler::optimizeloops))
        //.def("optimizeloops",&optimizeloops)
        .def("potenergy",&spl::LoopSampler::potenergy)
        .def("loopenes",&loopenes)
        .def("loopatomcrds",&loopatomcrds)
        .def("regenloops",&regenloops)
        //.def("regenuseloops",&regenuseloops)
        //.def("regenuseloops",[](spl::LoopSampler& tl) -> std::vector<NSPgeometry::XYZ>& {
        //    return *tl.regenuseloops();}, py::return_value_policy::reference_internal)
        .def("updatelooppools",&spl::LoopSampler::updatelooppools)
        .def("readparams",&spl::LoopSampler::readparams)
        .def("useparams",&spl::LoopSampler::useparams)
        .def("runMCCycles",&spl::LoopSampler::runMCCycles)
        .def_static("inloopblks",&spl::LoopSampler::inloopblks)
        .def("mypara", overload_cast_<>()(&spl::LoopSampler::mypara))
        .def("mypara", overload_cast_<>()(&spl::LoopSampler::mypara,py::const_))
        .def("imol", overload_cast_<spl::LoopSampler>()(&imol))
        .def("reindexedflankingsites", overload_cast_<>()(&spl::LoopSampler::reindexedflankingsites))
        .def("reindexedflankingsites", overload_cast_<>()(&spl::LoopSampler::reindexedflankingsites,py::const_));
    m.def("makelpsamplerparam",&spl::makelpsamplerparam);
    m.def("readpara_lpsampler",&spl::readpara_lpsampler);
    m.def("mkloopsampler",&spl::mkloopsampler);
    //molsdrun.h
    py::class_<spl::EmptySDCallBack>(m, "EmptySDCallBack")
        .def("__call__", &spl::EmptySDCallBack::operator());
    py::class_<spl::MolSDRun>(m, "MolSDRun")
        .def(py::init<>())
        .def("initrandomengine",&spl::MolSDRun::initrandomengine<unsigned int>)
        .def("sd",&sd)
        .def("imol",overload_cast_<spl::MolSDRun>()(&imol))
        .def("shakebds",&shakebds)
        // .def("init",&spl::MolSDRun::init)
        .def("init",[](MolSDRun &sd,SDRunPara para, NSPintrct::IntrctMol imol, bool shakeinit=true){
            std::shared_ptr<NSPintrct::IntrctMol> newimol=std::shared_ptr<NSPintrct::IntrctMol>(new IntrctMol(imol));
            sd.init(para,newimol,shakeinit);
        })
        .def("addrestraint",&spl::MolSDRun::addrestraint)
        .def("addrestraint_contact",[](MolSDRun &sd, ContactRestraint res){
            std::shared_ptr<NSPsampling::ResTraint> cr=std::shared_ptr<NSPsampling::ResTraint>(new NSPsampling::ContactRestraint(res));
            sd.addrestraint(cr);
        })
        .def("runsteps",&spl::MolSDRun::runsteps<spl::EmptySDCallBack>)
        .def("runsteps",&spl::MolSDRun::runsteps<spl::SDCallBack_default>)
        .def("state", overload_cast_<>()(&spl::MolSDRun::state))
        .def("state", overload_cast_<>()(&spl::MolSDRun::state,py::const_))
        .def("reinitstate",&spl::MolSDRun::reinitstate)
        .def("nstepsrun",&spl::MolSDRun::nstepsrun)
        .def("temperature",&spl::MolSDRun::temperature)
        .def("shakeon",overload_cast_<>()(&spl::MolSDRun::shakeon))
        .def("shakeon",overload_cast_<>()(&spl::MolSDRun::shakeon,py::const_))
        .def("nblsteps",overload_cast_<>()(&spl::MolSDRun::nblsteps))
        .def("nblsteps",overload_cast_<>()(&spl::MolSDRun::nblsteps,py::const_))
        .def("changetemperature", overload_cast_<double>()(&spl::MolSDRun::changetemperature))
        .def("changetemperature", overload_cast_<double,const std::vector<int> &>()(&spl::MolSDRun::changetemperature))
        .def("changetemperature", overload_cast_<double,int>()(&spl::MolSDRun::changetemperature))
        //.def("temperaturegroups_sp", overload_cast_<>()(&spl::MolSDRun::temperaturegroups_sp))
        //.def("temperaturegroups_sp", overload_cast_<>()(&spl::MolSDRun::temperaturegroups_sp,py::const_))
        .def("bathtemperatures",&spl::MolSDRun::bathtemperatures)
        .def("restraints",&spl::MolSDRun::restraints)
        .def("initpara",&spl::MolSDRun::initpara);
    m.def("newmolsdrun",&spl::newmolsdrun);
    m.def("makemolsdrun",&spl::makemolsdrun);
    m.def("restraints_energies",[](MolSDRun sdrun){
        std::vector<double> restraintenes;
        for(auto e:sdrun.restraints()) restraintenes.push_back(e->ene());
        return restraintenes;
    });


    //sdrunpara.h
    py::class_<spl::SDRunPara> (m, "SDRunPara")
        .def(py::init<>())
        .def(py::init<const std::vector<std::string> &>())
        .def_readwrite("jobname", &spl::SDRunPara::jobname)
        .def_readwrite("outpdbfile", &spl::SDRunPara::outpdbfile)
        .def_readwrite("topnpdbname", &spl::SDRunPara::topnpdbname)
        .def_readwrite("storetopn", &spl::SDRunPara::storetopn)
        .def_readwrite("emaxforstore", &spl::SDRunPara::emaxforstore)
        .def_readwrite("nr_rmsd", &spl::SDRunPara::nr_rmsd)
        .def_readwrite("epot_decay", &spl::SDRunPara::epot_decay)
        .def_readwrite("restraintsfile", &spl::SDRunPara::restraintsfile)
        .def_readwrite("doshake", &spl::SDRunPara::doshake)
        .def_readwrite("shakeinitcrd", &spl::SDRunPara::shakeinitcrd)
        .def_readwrite("newnbliststeps", &spl::SDRunPara::newnbliststeps)
        .def_readwrite("newsscodesteps", &spl::SDRunPara::newsscodesteps)
        .def_readwrite("printsteps", &spl::SDRunPara::printsteps)
        .def_readwrite("savepdbsteps", &spl::SDRunPara::savepdbsteps)
        .def_readwrite("gamma", &spl::SDRunPara::gamma)
        .def_readwrite("atommass", &spl::SDRunPara::atommass)
        .def_readwrite("timestep", &spl::SDRunPara::timestep)
        .def_readwrite("randomseed", &spl::SDRunPara::randomseed)
        .def_readwrite("doannealing", &spl::SDRunPara::doannealing)
        .def_readwrite("atomgroups", &spl::SDRunPara::atomgroups)
        .def_readwrite("temperaturegrps", &spl::SDRunPara::temperaturegrps)
        .def_readwrite("temperatures", &spl::SDRunPara::temperatures)
        .def_readwrite("annealingscheme", &spl::SDRunPara::annealingscheme)
        .def_readwrite("annealinggroup", &spl::SDRunPara::annealinggroup);
    m.def("makesdrunparam",&spl::makesdrunparam);
    //runsddefaults.h
    py::class_<spl::SDCallBack_default> sdcallback_default(m, "SDCallBack_default");
    sdcallback_default.def(py::init<const spl::SDRunPara &>())
        .def_readwrite("printsteps", &spl::SDCallBack_default::printsteps)
        .def_readwrite("savepdbsteps", &spl::SDCallBack_default::savepdbsteps)
        .def_readwrite("doannealing", &spl::SDCallBack_default::doannealing)
        .def_readwrite("pdbfilename", &spl::SDCallBack_default::pdbfilename)
        .def_readwrite("storedconfig", &spl::SDCallBack_default::storedconfig)
        .def_readwrite("epot_av", &spl::SDCallBack_default::epot_av)
        .def_readwrite("lastwritetopstep", &spl::SDCallBack_default::lastwritetopstep)
        .def("__call__", &spl::SDCallBack_default::operator());
    py::class_<spl::SDCallBack_default::Config> (sdcallback_default, "Config")
        .def_readwrite("crd", &spl::SDCallBack_default::Config::crd)
        .def_readwrite("energies", &spl::SDCallBack_default::Config::energies)
        .def_readwrite("restraint_enes", &spl::SDCallBack_default::Config::restraint_enes)
        .def_readwrite("pot_tot", &spl::SDCallBack_default::Config::pot_tot)
        .def_readwrite("cfstep", &spl::SDCallBack_default::Config::cfstep);
    m.def("sdreadparameters",&spl::readparameters);



    //pocketdocker.h  cyx
    py::class_<spl::NbrFndrGrid>(m, "NbrFndrGrid")
        .def_readwrite("grid", &spl::NbrFndrGrid::grid)
        .def_readwrite("gridneighbors", &spl::NbrFndrGrid::gridneighbors)
        .def_readwrite("cutoff2", &spl::NbrFndrGrid::cutoff2)
        .def("setgrid", &spl::NbrFndrGrid::setgrid)
        .def("setgridneighbors", &spl::NbrFndrGrid::setgridneighbors)
        .def("candidateneighbors", &spl::NbrFndrGrid::candidateneighbors);
    py::class_<spl::BckBnStMatcher>(m, "BckBnStMatcher")
        .def_static("sitedistance2", &spl::BckBnStMatcher::sitedistance2)
        .def("matchsite", &spl::BckBnStMatcher::matchsite)
        .def_readwrite("candidatesites", &spl::BckBnStMatcher::candidatesites)
        .def_readwrite("maxdist2", &spl::BckBnStMatcher::maxdist2);
    py::class_<spl::PocketMover>(m, "PocketMover")
        .def_readwrite("ligandcentersphere", &spl::PocketMover::ligandcentersphere)
        .def_readwrite("maxstep_trans", &spl::PocketMover::maxstep_trans)
        .def_readwrite("maxstep_rotate", &spl::PocketMover::maxstep_rotate)
        .def_readwrite("prigidmv", &spl::PocketMover::prigidmv)
        .def_readwrite("refcrd", &spl::PocketMover::refcrd)
        .def_readwrite("crd_before", &spl::PocketMover::crd_before)
        .def_readwrite("movencomponents", &spl::PocketMover::movencomponents)
        .def_readwrite("pocketclashes", &spl::PocketMover::pocketclashes)
        .def("stepback", &spl::PocketMover::stepback)
        .def("randommv", &spl::PocketMover::randommv)
        .def("mvtosphere", &spl::PocketMover::mvtosphere)
        .def("mvtofit", &spl::PocketMover::mvtofit);
    py::class_<spl::PckDckrParam> pckdckrparam(m, "PckDckrParam");
    pckdckrparam.def(py::init<>())
        .def_readwrite("proteinpdb", &spl::PckDckrParam::proteinpdb)
        .def_readwrite("pocketfile", &spl::PckDckrParam::pocketfile)
        .def_readwrite("ligandlocationmode", &spl::PckDckrParam::ligandlocationmode)
        .def_readwrite("ligandsphere_center", &spl::PckDckrParam::ligandsphere_center)
        .def_readwrite("ligandcenteratoms", &spl::PckDckrParam::ligandcenteratoms)
        .def_readwrite("ligandsphereradius", &spl::PckDckrParam::ligandsphereradius)
        .def_readwrite("gridcenter", &spl::PckDckrParam::gridcenter)
        .def_readwrite("gridlengths", &spl::PckDckrParam::gridlengths)
        .def_readwrite("gridstep", &spl::PckDckrParam::gridstep)
        .def_readwrite("maxtranslation", &spl::PckDckrParam::maxtranslation)
        .def_readwrite("maxrotation", &spl::PckDckrParam::maxrotation);
    py::enum_<spl::PckDckrParam::POCKETLOCATIONMODE>(pckdckrparam, "POCKETLOCATIONMODE")
        .value("CENTERXYZ", spl::PckDckrParam::POCKETLOCATIONMODE::CENTERXYZ)
        .value("CENTERATOMS", spl::PckDckrParam::POCKETLOCATIONMODE::CENTERATOMS)
        .value("REFATOMS", spl::PckDckrParam::POCKETLOCATIONMODE::REFATOMS)
        .export_values();
    //m.def("mk_pckdckrparam",&spl::mk_pckdckrparam);  
    py::class_<spl::PocketDocker>(m, "PocketDocker")
        //.def("setup",&spl::PocketDocker::setup)
        .def("micromove", &spl::PocketDocker::micromove);
    //screen.h
    py::class_<Screen::SpherePoints>(m, "SpherePoints")
        .def(py::init<>())
        .def_readwrite("points", &Screen::SpherePoints::points);
    m.def("expose", &Screen::expose);
    m.def("formhb", &Screen::formhb);
    //savedloops.h
    typedef spl::LoopSampler::FlankingSite FlankingSite;
    py::class_<spl::SavedLoops> savedloops(m, "SavedLoops");
    savedloops.def(py::init<>())
        .def_readwrite("nloops", &spl::SavedLoops::nloops)
        .def_readwrite("nloopatoms", &spl::SavedLoops::nloopatoms)
        .def_readwrite("min_nloopatom", &spl::SavedLoops::min_nloopatom)
        .def_readwrite("reindexedflankingsites", &spl::SavedLoops::reindexedflankingsites)
        .def_readwrite("looplengths", &spl::SavedLoops::looplengths)
        .def_readwrite("loopsequences", &spl::SavedLoops::loopsequences)
        .def_readwrite("loopfile", &spl::SavedLoops::loopfile)
        .def_readwrite("imol", &spl::SavedLoops::imol)
        .def_readwrite("topnconfigs", &spl::SavedLoops::topnconfigs)
        .def_readwrite("configs", &spl::SavedLoops::configs)
        .def_readwrite("startconfig", &spl::SavedLoops::startconfig)
        .def_readwrite("pdbconfig", &spl::SavedLoops::pdbconfig)
        .def_static("rmsd", &spl::SavedLoops::rmsd)
        .def_static("pprmsd", &spl::SavedLoops::pprmsd)
        .def("rmsdcut", &spl::SavedLoops::rmsdcut)
        .def("rmsdcut_de", &spl::SavedLoops::rmsdcut_de)
        .def("init", &spl::SavedLoops::init)
        .def("readnextconfig", &spl::SavedLoops::readnextconfig)
        .def("copypdbconfig", &spl::SavedLoops::copypdbconfig)
        .def("readstartconfig", &spl::SavedLoops::readstartconfig)
        .def("readallconfigs", &spl::SavedLoops::readallconfigs)
        .def("determinetopn", &spl::SavedLoops::determinetopn)
        .def("determinetopnPP", &spl::SavedLoops::determinetopnPP)
        .def("savenonredundant", &spl::SavedLoops::savenonredundant)
        .def("writepdb", &spl::SavedLoops::writepdb)
        .def("writepdb_topn", &spl::SavedLoops::writepdb_topn)
        .def("writermsdene", &spl::SavedLoops::writermsdene);
    py::class_<spl::SavedLoops::Config>(savedloops, "Config")
        .def_readwrite("energies", &spl::SavedLoops::Config::energies)
        .def_readwrite("crds", &spl::SavedLoops::Config::crds)
        .def_readwrite("id", &spl::SavedLoops::Config::id);
    py::class_<spl::SavedLoops::CompareConfigPP>(savedloops, "CompareConfigPP")
        .def_readwrite("slp", &spl::SavedLoops::CompareConfigPP::slp)
        .def("__call__", &spl::SavedLoops::CompareConfigPP::operator());
    py::class_<spl::SavedLoops::CompareConfig>(savedloops, "CompareConfig")
        .def_readwrite("slp", &spl::SavedLoops::CompareConfig::slp)
        .def("__call__", &spl::SavedLoops::CompareConfig::operator());
}
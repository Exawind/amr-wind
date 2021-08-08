#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;

static
void 
print_usage (int,
             char* argv[])
{
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " infile=f1 [options] \n\tOptions:\n";
  exit(1);
}

std::string
getFileRoot(const std::string& infile)
{
  const auto& tokens = Tokenize(infile,std::string("/"));
  return tokens[tokens.size()-1];
}

int
main (int   argc,
      char* argv[])
{
  Initialize(argc,argv);
  {
    if (argc < 2)
      print_usage(argc,argv);

    ParmParse pp;

    if (pp.contains("help"))
      print_usage(argc,argv);

    // Open plotfile header and create an amrData object pointing into it
    std::string plotFileName; pp.get("infile",plotFileName);

    PlotFileData pf(plotFileName);
    int finestLevel = pf.finestLevel();
    int Nlev = 1; pp.query("levels",Nlev);
    int ratio = 2; pp.query("ratio",ratio);

    int nComp, sComp;
    const auto& varNames = pf.varNames();
    Vector<std::string> varnames(1);
    if (pp.countval("varname") != 0) {
      std::string varname;
      pp.get("varname",varname);
      nComp = 1;
      sComp = -1;
      for (int j=0; j<varNames.size() && sComp==-1; ++j) {
        if (varNames[j] == varname) sComp=j;
      }
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(sComp!=-1,"variable not found in plotfile");
      varnames[0] = varname;
    } else {
      nComp = varNames.size();
      sComp = 0;
      varnames = varNames;
    }

    Vector<MultiFab *> data(Nlev);
    int coordSys = pf.coordSys();
    const auto problo = pf.probLo();
    const auto probhi = pf.probHi();
    RealBox rb(problo,probhi);
    Array<int,AMREX_SPACEDIM> is_per = {{0,0,1}};
    Vector<int> levelSteps(Nlev);
    Vector<IntVect> refRatio(Nlev-1);
    Vector<Geometry> geoms(Nlev);
    Box fpd = pf.probDomain(0);
    Geometry fgeom(fpd,rb,coordSys,is_per);
    auto mf = (nComp == 1 ? pf.get(0,varnames[0]) : pf.get(0));
    const auto& fba = mf.boxArray();
    const auto& dm = mf.DistributionMap();

    for (int lev=0; lev<Nlev; ++lev) {

	  int r = std::pow(ratio,Nlev-lev);
	  std::cout << "level " << lev << "; ratio from original = " << r << std::endl;

      Box cpd = Box(fpd).coarsen(r);
      geoms[lev].define(cpd,rb,coordSys,is_per);
      levelSteps[lev] = pf.levelStep(0);
      if (lev < Nlev-1) {
        refRatio[lev] = IntVect(AMREX_D_DECL(ratio,ratio,ratio));
      }

      BoxArray cba = BoxArray(fba).coarsen(r);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(BoxArray(cba).refine(r) == fba,"BoxArray doesn't coarsen correctly");

      data[lev] = new MultiFab(cba,dm,nComp,0);
      average_down(mf,*data[lev],fgeom,geoms[lev],0,nComp,r);

    }

    Vector<const MultiFab*> dat(Nlev);
    for (int i=0; i<Nlev; ++i) {
      dat[i] = data[i];
    }
    WriteMultiLevelPlotfile(getFileRoot(plotFileName)+"_crse",Nlev,dat,varnames,geoms,pf.time(),levelSteps,refRatio);
  }
  Finalize();
  return 0;
}

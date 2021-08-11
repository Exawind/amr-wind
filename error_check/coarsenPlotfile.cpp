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
        std::string plotFileName1; pp.get("referencefile",plotFileName1);
        std::string plotFileName2; pp.get("amrfile",plotFileName2);
  
        PlotFileData pf1(plotFileName1);
        PlotFileData pf2(plotFileName2);
        int finest_level = pf2.finestLevel();
        int Nlev = finest_level + 1;
        //int Nlev = 1; pp.query("levels",Nlev);
        int ratio = 2;// pp.query("ratio",ratio);
        int referenceRefinement = 2;// pp.query("ratio",ratio);
  
        const auto& dx_1 = pf1.cellSize(0);
        const auto& dx_2 = pf2.cellSize(Nlev-1);
  	    std::cout << "dx_1 = " << dx_1 << "; dx_2 = " << dx_2 << std::endl;
        bool not_match = AMREX_D_TERM(   dx_1[0] * referenceRefinement != dx_2[0],
                                      || dx_1[1] * referenceRefinement != dx_2[1],
                                      || dx_1[2] * referenceRefinement != dx_2[2] );
        if (not_match) {
            amrex::Abort("ERROR: grid dx at finest level does not match");
        }
  
        int nComp, sComp;
        const auto& varNames = pf1.varNames();
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
  
  	    for (int j=0; j<varnames.size(); ++j){
  	    	std::cout << "varnames[" << j << "] = " << varnames[j] << std::endl;
  	    }
  
        Vector<MultiFab *> data(Nlev);
        int coordSys = pf1.coordSys();
        const auto problo = pf1.probLo();
        const auto probhi = pf1.probHi();
        RealBox rb(problo,probhi);
        Array<int,AMREX_SPACEDIM> is_per = {{0,0,1}};
        Vector<int> levelSteps(Nlev);
        Vector<IntVect> refRatio(Nlev-1);
        Vector<Geometry> geoms(Nlev);
        Box fpd = pf1.probDomain(0);
        Geometry fgeom(fpd,rb,coordSys,is_per);
        auto mf1 = (nComp == 1 ? pf1.get(0,varnames[0]) : pf1.get(0));
        const auto& fba = mf1.boxArray();
        const auto& dm = mf1.DistributionMap();
  
  
  	    amrex::Real L2error = 0.0;
        for (int lev=0; lev<=finest_level; ++lev) {
  
  	        int r = std::pow(ratio,finest_level-lev) * referenceRefinement;
  	        std::cout << "coarsening reference solution (level " << Nlev  << ") to level " << lev << "; ratio from reference = " << r << std::endl;
  
            Box cpd = Box(fpd).coarsen(r);
            geoms[lev].define(cpd,rb,coordSys,is_per);
            levelSteps[lev] = pf1.levelStep(0);
            refRatio[lev] = IntVect(AMREX_D_DECL(ratio,ratio,ratio));
  
            BoxArray cba = BoxArray(fba).coarsen(r);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(BoxArray(cba).refine(r) == fba,"BoxArray doesn't coarsen correctly");
  
            data[lev] = new MultiFab(cba,dm,nComp,0);
            average_down(mf1,*data[lev],fgeom,geoms[lev],0,nComp,r);
            
			auto mf2_lev = (nComp == 1 ? pf2.get(lev,varnames[0]) : pf2.get(lev));

			MultiFab mf_b;
		    mf_b.define(mf2_lev.boxArray(), mf2_lev.DistributionMap(), 1, 0);
			mf_b.ParallelCopy(*data[lev]);

			amrex::iMultiFab level_mask;
       	    if (lev < finest_level) {
				auto mf2_levp1 = (nComp == 1 ? pf2.get(lev + 1,varnames[0]) : pf2.get(lev + 1));
       	        level_mask = makeFineMask(
       	            mf2_lev.boxArray(), mf2_lev.DistributionMap(),
       	            mf2_levp1.boxArray(), amrex::IntVect(2), 1, 0);
       	    } else {
       	        level_mask.define(
       	            mf2_lev.boxArray(), mf2_lev.DistributionMap(),
       	            1, 0, amrex::MFInfo());
       	        level_mask.setVal(1);
       	    }	

  	    	// compare with refined solution
  	    	
            const auto& dx = pf2.cellSize(lev);
  	        for (amrex::MFIter mfi(mf2_lev); mfi.isValid(); ++mfi) {
                const auto& vbx = mfi.validbox();
  	            std::cout << "box: " << vbx << std::endl;
  	            auto mask = level_mask.array(mfi); 
			    auto var1 = mf_b.array(mfi);
  	            auto var2 = mf2_lev.array(mfi);
  	            amrex::Loop(
                      vbx, [=, &L2error] (int i, int j, int k) noexcept {
						  const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                		  const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                		  const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

                          const amrex::Real diff = mask(i, j ,k) * (var2(i, j, k, 0) - var1(i, j, k, 0));
                          L2error += std::pow(diff,2);

						 // if(diff != 0.0){
						 // 	  std::cout << "y + z = " << y + z << std::endl;
						 //     std::cout << "var1 = " << var1(i,j,k,0) << "; var2 = " << var2(i,j,k,0) << std::endl;
						 // }
                      });
  	        }
  	        std::cout << "level " << lev << " done; L2error squared = " << L2error << std::endl;
  
        }

		L2error = std::sqrt(L2error);
  		std::cout << "L2error = " << L2error << std::endl;	
  
        Vector<const MultiFab*> dat(Nlev);
        for (int i=0; i<Nlev; ++i) {
          dat[i] = data[i];
        }
        WriteMultiLevelPlotfile(getFileRoot(plotFileName1)+"_crse",Nlev,dat,varnames,geoms,pf1.time(),levelSteps,refRatio);
    }
    Finalize();
    return 0;
}

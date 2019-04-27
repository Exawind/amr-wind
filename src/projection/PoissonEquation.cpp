#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>

#include <PoissonEquation.H>
#include <projection_F.H>

using namespace amrex;

//
// Constructor: 
// We set up everything which doesn't change between timesteps here
//
PoissonEquation::PoissonEquation(AmrCore* _amrcore,
                                 Vector<std::unique_ptr<EBFArrayBoxFactory>>* _ebfactory,
                                 Vector<std::unique_ptr<IArrayBox>>& bc_ilo, 
                                 Vector<std::unique_ptr<IArrayBox>>& bc_ihi, 
                                 Vector<std::unique_ptr<IArrayBox>>& bc_jlo, 
                                 Vector<std::unique_ptr<IArrayBox>>& bc_jhi, 
                                 Vector<std::unique_ptr<IArrayBox>>& bc_klo, 
                                 Vector<std::unique_ptr<IArrayBox>>& bc_khi,
                                 int _nghost)
{
    // Get inputs from ParmParse
	readParameters();
    
    if(verbose > 0)
    {
        amrex::Print() << "Constructing PoissonEquation class" << std::endl;
    }

    // Set AmrCore and ebfactory based on input, fetch some data needed in constructor
    amrcore = _amrcore; 
    ebfactory = _ebfactory;
    nghost = _nghost;
    Vector<Geometry> geom = amrcore->Geom();
    Vector<BoxArray> grids = amrcore->boxArray();
    Vector<DistributionMapping> dmap = amrcore->DistributionMap();
    int max_level = amrcore->maxLevel();
    
    // Whole domain
    Box domain(geom[0].Domain());

    // The boundary conditions need only be set at level 0
    set_ppe_bc(bc_lo, bc_hi,
               domain.loVect(), domain.hiVect(),
			   &nghost,
               bc_ilo[0]->dataPtr(), bc_ihi[0]->dataPtr(),
               bc_jlo[0]->dataPtr(), bc_jhi[0]->dataPtr(),
               bc_klo[0]->dataPtr(), bc_khi[0]->dataPtr());

    // Resize and reset sigma
    sigma.resize(max_level + 1);
    for(int lev = 0; lev <= max_level; lev++)
    {
        sigma[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost, 
                                      MFInfo(), *(*ebfactory)[lev]));
    }

	// First define the matrix.
    // Class MLNodeLaplacian describes the following operator:
    //
    //       del dot (sigma grad) phi = rhs,
    //
    // where phi and rhs are nodal, and sigma is cell-centered
	LPInfo info;
	info.setMaxCoarseningLevel(mg_max_coarsening_level);

    matrix.define(geom, grids, dmap, info, GetVecOfConstPtrs(*ebfactory));

    matrix.setGaussSeidel(true);
    matrix.setHarmonicAverage(false);

	// LinOpBCType Definitions are in amrex/Src/Boundary/AMReX_LO_BCTYPES.H
	matrix.setDomainBC
    (
        {(LinOpBCType) bc_lo[0], (LinOpBCType) bc_lo[1], (LinOpBCType) bc_lo[2]},
        {(LinOpBCType) bc_hi[0], (LinOpBCType) bc_hi[1], (LinOpBCType) bc_hi[2]}
    );
}

PoissonEquation::~PoissonEquation()
{
}

void PoissonEquation::readParameters()
{
    ParmParse pp("projection");

    pp.query("verbose", verbose);
    pp.query("mg_verbose", mg_verbose);
    pp.query("mg_cg_verbose", mg_cg_verbose);
    pp.query("mg_max_iter", mg_max_iter);
    pp.query("mg_cg_maxiter", mg_cg_maxiter);
    pp.query("mg_max_fmg_iter", mg_max_fmg_iter);
    pp.query("mg_max_coarsening_level", mg_max_coarsening_level);
    pp.query("mg_rtol", mg_rtol);
    pp.query("mg_atol", mg_atol);
    pp.query( "bottom_solver_type", bottom_solver_type);
}

void PoissonEquation::updateInternals(AmrCore* amrcore_in, 
                                      Vector<std::unique_ptr<EBFArrayBoxFactory>>* ebfactory_in)
{
    // This must be implemented when we want dynamic meshing
    //
    amrex::Print() << "ERROR: PoissonEquation::updateInternals() not yet implemented" << std::endl;
    amrex::Abort(); 
}

// 
// Set the user-supplied settings for the MLMG solver
// (this must be done every time step, since MLMG is created after updating matrix
//
void PoissonEquation::setSolverSettings(MLMG& solver)
{
    // The default bottom solver is now bicgcg
    if (bottom_solver_type == "smoother")
    {
       solver.setBottomSolver(MLMG::BottomSolver::smoother);
    }
    else if (bottom_solver_type == "bicg")
    {
       solver.setBottomSolver(MLMG::BottomSolver::bicgstab);
    }
    else if (bottom_solver_type == "cg")
    {
       solver.setBottomSolver(MLMG::BottomSolver::cg);
    }
    else if (bottom_solver_type == "bicgcg")
    {
       solver.setBottomSolver(MLMG::BottomSolver::bicgcg);
    }
    else if (bottom_solver_type == "cgbicg")
    {
       solver.setBottomSolver(MLMG::BottomSolver::cgbicg);
    }
    else if (bottom_solver_type == "hypre")
    {
       solver.setBottomSolver(MLMG::BottomSolver::hypre);
    }

    // Maximum iterations for MultiGrid / ConjugateGradients
	solver.setMaxIter(mg_max_iter);
	solver.setMaxFmgIter(mg_max_fmg_iter);
	solver.setCGMaxIter(mg_cg_maxiter);

    // Verbosity for MultiGrid / ConjugateGradients
	solver.setVerbose(mg_verbose);
	solver.setCGVerbose(mg_cg_verbose);

	// This ensures that ghost cells of phi are correctly filled when returned from the solver
	solver.setFinalFillBC(true);
}

//
// Solve Poisson Equation:
//
//                  div( 1/rho * grad(phi) ) = div(u)
//
// We output grad(phi) / rho into "fluxes"
//
void PoissonEquation::solve(Vector<std::unique_ptr<MultiFab>>& phi,
			                Vector<std::unique_ptr<MultiFab>>& fluxes,
                            const Vector<std::unique_ptr<MultiFab>>& ro, 
                            const Vector<std::unique_ptr<MultiFab>>& divu)
{
    for(int lev = 0; lev <= amrcore->finestLevel(); lev++)
    {
        // Set the coefficients to equal 1 / ro 
        sigma[lev]->setVal(1.0);
        MultiFab::Divide(*sigma[lev], *ro[lev], 0, 0, 1, nghost);
        matrix.setSigma(lev, *sigma[lev]);

        // By this point we must have filled the Dirichlet values of phi in ghost cells
        matrix.setLevelBC(lev, GetVecOfConstPtrs(phi)[lev]);
    }

    // Set up the solver
	MLMG solver(matrix);
    setSolverSettings(solver);

    // Solve!
	solver.solve(GetVecOfPtrs(phi), GetVecOfConstPtrs(divu), mg_rtol, mg_atol);

    // Get fluxes (grad(phi) / rho)
    solver.getFluxes(amrex::GetVecOfPtrs(fluxes));
}


#include <AMReX_MultiFabUtil.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EB_utils.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>

#include <DiffusionOp.H>
#include <diffusion_F.H>

using namespace amrex;

// Define unit vectors to easily convert indices
extern const amrex::IntVect e_x;
extern const amrex::IntVect e_y;
extern const amrex::IntVect e_z;

//
// Constructor:
// We set up everything which doesn't change between timesteps here
//
DiffusionOp::DiffusionOp(AmrCore* _amrcore,
                         Vector<std::unique_ptr<EBFArrayBoxFactory>>* _ebfactory,
                         std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_bc_lo,
                         std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_bc_hi,
                         int _nghost)
{
    if(verbose > 0)
        amrex::Print() << "Constructing DiffusionOp class" << std::endl;

    nghost = _nghost;

    m_bc_lo = a_bc_lo;
    m_bc_hi = a_bc_hi;

    // Get inputs from ParmParse
    readParameters();

    // Actually do the setup work here
    setup(_amrcore, _ebfactory);
}

void DiffusionOp::setup(AmrCore* _amrcore, 
                        Vector<std::unique_ptr<EBFArrayBoxFactory>>* _ebfactory)
{
    // The amrcore boxArray and DistributionMap change when we regrid so we must pass the new object in here.
    amrcore = _amrcore;

    // The ebfactory changes when we regrid so we must pass it in here.
    ebfactory = _ebfactory;

    geom  = amrcore->Geom();
    grids = amrcore->boxArray();
    dmap  = amrcore->DistributionMap();

    max_level = amrcore->maxLevel();

    // Resize and reset data
    b.resize(max_level + 1);
    phi.resize(max_level + 1);
    rhs.resize(max_level + 1);
    vel_eb.resize(max_level + 1);

    for(int lev = 0; lev <= max_level; lev++)
    {
        for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
        {
            BoxArray edge_ba = grids[lev];
            edge_ba.surroundingNodes(dir);
            b[lev][dir].reset(new MultiFab(edge_ba, dmap[lev], 1, nghost,
                                           MFInfo(), *(*ebfactory)[lev]));
        }
        phi[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, 1,
                                    MFInfo(), *(*ebfactory)[lev]));

        // No ghost cells needed for rhs
        rhs[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, 0,
                                    MFInfo(), *(*ebfactory)[lev]));

        vel_eb[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, nghost,
                                       MFInfo(), *(*ebfactory)[lev]));
        vel_eb[lev]->setVal(0.0);
    }

    // Define the matrix.
    LPInfo info;
    info.setMaxCoarseningLevel(mg_max_coarsening_level);
    matrix.reset(new MLEBTensorOp(geom, grids, dmap, info, GetVecOfConstPtrs(*ebfactory)));

    // It is essential that we set MaxOrder to 2 if we want to use the standard
    // phi(i)-phi(i-1) approximation for the gradient at Dirichlet boundaries.
    // The solver's default order is 3 and this uses three points for the gradient.
    matrix->setMaxOrder(2);

    // LinOpBCType Definitions are in amrex/Src/Boundary/AMReX_LO_BCTYPES.H
    matrix->setDomainBC(m_bc_lo, m_bc_hi);
}

DiffusionOp::~DiffusionOp()
{
}

void DiffusionOp::readParameters()
{
    ParmParse pp("diffusion");

    pp.query("verbose", verbose);
    pp.query("mg_verbose", mg_verbose);
    pp.query("mg_cg_verbose", mg_cg_verbose);
    pp.query("mg_max_iter", mg_max_iter);
    pp.query("mg_cg_maxiter", mg_cg_maxiter);
    pp.query("mg_max_fmg_iter", mg_max_fmg_iter);
    pp.query("mg_max_coarsening_level", mg_max_coarsening_level);
    pp.query("mg_rtol", mg_rtol);
    pp.query("mg_atol", mg_atol);
    pp.query("bottom_solver_type", bottom_solver_type);
}

//
// Solve the matrix equation
//
void DiffusionOp::solve(      Vector<std::unique_ptr<MultiFab>>& vel_in,
                        const Vector<std::unique_ptr<MultiFab>>& ro_in,
                        const Vector<std::unique_ptr<MultiFab>>& mu_in,
                        Real dt)
{
    BL_PROFILE("DiffusionOp::solve");

    // Update the coefficients of the matrix going into the solve based on the current state of the
    // simulation. Recall that the relevant matrix is
    //
    //      alpha a - beta div ( b grad )   <--->   rho - dt div ( mu grad )
    //
    // So the constants and variable coefficients are:
    //
    //      alpha: 1
    //      beta: dt
    //      a: ro
    //      b: mu

    // Set alpha and beta
    matrix->setScalars(1.0, dt);

    for(int lev = 0; lev <= max_level; lev++)
    {
        // Compute the spatially varying b coefficients (on faces) to equal the apparent viscosity
        average_cellcenter_to_face(GetArrOfPtrs(b[lev]), *mu_in[lev], geom[lev]);
        for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
            b[lev][dir]->FillBoundary(geom[lev].periodicity());
        
        // This sets the coefficients
        matrix->setACoeffs(lev, (*ro_in[lev]));
        matrix->setShearViscosity  (lev, GetArrOfConstPtrs(b[lev]));
        matrix->setEBShearViscosity(lev, (*mu_in[lev]));
    }

    if(verbose > 0)
        amrex::Print() << "Diffusing velocity components all together..." << std::endl; 

    for(int lev = 0; lev <= max_level; lev++)
    {
        // Set the right hand side to equal rho
        MultiFab::Copy((*rhs[lev]),(*vel_in[lev]), 0, 0, AMREX_SPACEDIM, 0);

        // Multiply rhs by rho to get momentum
        // Note that vel holds the updated velocity:
        //
        //      u_old + dt ( - u grad u + div ( mu (grad u)^T ) / rho - grad p / rho + gravity )
        //
        for (int i = 0; i < 3; i++)
           MultiFab::Multiply((*rhs[lev]), (*ro_in[lev]), 0, i, 1, 0);

        // By this point we must have filled the Dirichlet values of phi stored in ghost cells
        MultiFab::Copy(*phi[lev],*vel_in[lev], 0, 0, AMREX_SPACEDIM, 1);
        phi[lev]->FillBoundary(geom[lev].periodicity());
        matrix->setLevelBC(lev, GetVecOfConstPtrs(phi)[lev]);

        // matrix->setEBHomogDirichlet(lev, *mu_in[lev]);
    }

    MLMG solver(*matrix);
    setSolverSettings(solver);

    // This ensures that ghost cells of sol are correctly filled when returned from the solver
    solver.setFinalFillBC(true);

    solver.solve(GetVecOfPtrs(phi), GetVecOfConstPtrs(rhs), mg_rtol, mg_atol);

    for(int lev = 0; lev <= max_level; lev++)
    {
        phi[lev]->FillBoundary(geom[lev].periodicity());
        MultiFab::Copy(*vel_in[lev], *phi[lev], 0, 0, AMREX_SPACEDIM, 1);
    }

    if(verbose > 0)
        amrex::Print() << " done!" << std::endl;
}

//
// Set the user-supplied settings for the MLMG solver
// (this must be done every time step, since MLMG is created after updating matrix
//
void DiffusionOp::setSolverSettings(MLMG& solver)
{
    // The default bottom solver is BiCG
    if(bottom_solver_type == "smoother")
    {
        solver.setBottomSolver(MLMG::BottomSolver::smoother);
    }
    else if(bottom_solver_type == "hypre")
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

void DiffusionOp::ComputeDivTau(Vector<std::unique_ptr<MultiFab>>& divtau_out,
                                const Vector<std::unique_ptr<MultiFab>>& vel_in,
                                const Vector<std::unique_ptr<MultiFab>>& ro_in,
                                const Vector<std::unique_ptr<MultiFab>>& mu_in)
{
    BL_PROFILE("DiffusionOp::ComputeDivTau");

    Vector<std::unique_ptr<MultiFab> > divtau_aux(max_level+1);
    for(int lev = 0; lev <= max_level; lev++)
    {
       divtau_aux[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, nghost,
                                          MFInfo(), *(*ebfactory)[lev]));
       divtau_aux[lev]->setVal(0.0);
    }
 
    // Whole domain
    Box domain(geom[0].Domain());
 
    // We want to return div (mu grad)) phi
    matrix->setScalars(0.0, -1.0);
 
    // Compute the coefficients
    for (int lev = 0; lev <= max_level; lev++)
    {
        average_cellcenter_to_face( GetArrOfPtrs(b[lev]), *mu_in[lev], geom[lev] );
 
        for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
             b[lev][dir]->FillBoundary(geom[lev].periodicity());
 
        matrix->setShearViscosity  ( lev, GetArrOfConstPtrs(b[lev]));
        matrix->setEBShearViscosity( lev, (*mu_in[lev]));
        matrix->setLevelBC         ( lev, GetVecOfConstPtrs(vel_in)[lev] );
    }
 
    MLMG solver(*matrix);
 
    solver.apply(GetVecOfPtrs(divtau_aux), GetVecOfPtrs(vel_in));
 
    for(int lev = 0; lev <= max_level; lev++)
    {
       amrex::single_level_redistribute( lev, *divtau_aux[lev], *divtau_out[lev], 0, AMREX_SPACEDIM, geom);
 
       // Divide by density
       for (int n = 0; n < 3; n++)
           MultiFab::Divide( *divtau_out[lev], *ro_in[lev], 0, n, 1, 0 );
    }
}

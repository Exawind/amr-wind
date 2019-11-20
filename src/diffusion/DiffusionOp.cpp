#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EB_utils.H>
#endif

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
#ifdef AMREX_USE_EB
                         Vector<std::unique_ptr<EBFArrayBoxFactory>>* _ebfactory,
#endif
                         std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_velbc_lo,
                         std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_velbc_hi,
                         std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_scalbc_lo,
                         std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_scalbc_hi,
                         int _nghost,
                         int _probtype)
{
    if(verbose > 0)
        amrex::Print() << "Constructing DiffusionOp class" << std::endl;
    
    probtype = _probtype;
    nghost = _nghost;

    m_velbc_lo = a_velbc_lo;
    m_velbc_hi = a_velbc_hi;

    m_scalbc_lo = a_scalbc_lo;
    m_scalbc_hi = a_scalbc_hi;

    // Get inputs from ParmParse
    readParameters();

    // Actually do the setup work here
#ifdef AMREX_USE_EB
    setup(_amrcore, _ebfactory);
#else
    setup(_amrcore);
#endif
}

#ifdef AMREX_USE_EB
void DiffusionOp::setup(AmrCore* _amrcore, 
                        Vector<std::unique_ptr<EBFArrayBoxFactory>>* _ebfactory)
#else
void DiffusionOp::setup(AmrCore* _amrcore) 
#endif
{
    // The amrcore boxArray and DistributionMap change when we regrid so we must pass the new object in here.
    amrcore = _amrcore;

#ifdef AMREX_USE_EB
    // The ebfactory changes when we regrid so we must pass it in here.
    ebfactory = _ebfactory;
#endif

    geom  = amrcore->Geom();
    grids = amrcore->boxArray();
    dmap  = amrcore->DistributionMap();

    int max_level = amrcore->maxLevel();

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
#ifdef AMREX_USE_EB
            b[lev][dir].reset(new MultiFab(edge_ba, dmap[lev], 1, nghost,
                                           MFInfo(), *(*ebfactory)[lev]));
#else
            b[lev][dir].reset(new MultiFab(edge_ba, dmap[lev], 1, nghost, MFInfo()));
#endif
        }

#ifdef AMREX_USE_EB
        phi[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, 1,
                                    MFInfo(), *(*ebfactory)[lev]));
        rhs[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, 0,
                                    MFInfo(), *(*ebfactory)[lev]));
        vel_eb[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, nghost,
                                       MFInfo(), *(*ebfactory)[lev]));
#else
        phi[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, 1, MFInfo()));
        rhs[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo()));
        vel_eb[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, nghost, MFInfo()));
#endif
        vel_eb[lev]->setVal(0.0);
    }

    LPInfo info;
    info.setMaxCoarseningLevel(mg_max_coarsening_level);

    // 
    // Define the matrix for the viscous tensor solve.
    // 
#ifdef AMREX_USE_EB
     vel_matrix.reset(new MLEBTensorOp(geom, grids, dmap, info, GetVecOfConstPtrs(*ebfactory)));
#else
     vel_matrix.reset(new MLTensorOp(geom, grids, dmap, info));
#endif

    // It is essential that we set MaxOrder to 2 if we want to use the standard
    // phi(i)-phi(i-1) approximation for the gradient at Dirichlet boundaries.
    // The solver's default order is 3 and this uses three points for the gradient.
    vel_matrix->setMaxOrder(2);

    // LinOpBCType Definitions are in amrex/Src/Boundary/AMReX_LO_BCTYPES.H
    vel_matrix->setDomainBC(m_velbc_lo, m_velbc_hi);
    
    if(probtype == 35){
        amrex::Print() << " warning vel bcs are manually entered for probtype 35 in DiffusionOp.cpp " << std::endl;
        //fixme hard coded in
        Vector<Array<amrex::LinOpBCType,AMREX_SPACEDIM>> velbc_lo;
        Vector<Array<amrex::LinOpBCType,AMREX_SPACEDIM>> velbc_hi;

        velbc_lo.push_back({LinOpBCType::Periodic, LinOpBCType::Periodic, LinOpBCType::inhomogNeumann});
        velbc_lo.push_back({LinOpBCType::Periodic, LinOpBCType::Periodic, LinOpBCType::inhomogNeumann});
        velbc_lo.push_back({LinOpBCType::Periodic, LinOpBCType::Periodic, LinOpBCType::Dirichlet});

        velbc_hi.push_back({LinOpBCType::Periodic, LinOpBCType::Periodic, LinOpBCType::Dirichlet});
        velbc_hi.push_back({LinOpBCType::Periodic, LinOpBCType::Periodic, LinOpBCType::Dirichlet});
        velbc_hi.push_back({LinOpBCType::Periodic, LinOpBCType::Periodic, LinOpBCType::Dirichlet});

        vel_matrix->setDomainBC(velbc_lo,velbc_hi);
    }
    
    // 
    // Define the matrix for the scalar diffusion solve.
    // 
#ifdef AMREX_USE_EB
    scal_matrix.reset(new MLEBABecLap(geom, grids, dmap, info, GetVecOfConstPtrs(*ebfactory)));
#else
    scal_matrix.reset(new MLABecLaplacian(geom, grids, dmap, info));
#endif

    // It is essential that we set MaxOrder to 2 if we want to use the standard
    // phi(i)-phi(i-1) approximation for the gradient at Dirichlet boundaries.
    // The solver's default order is 3 and this uses three points for the gradient.
    scal_matrix->setMaxOrder(2);

    // LinOpBCType Definitions are in amrex/Src/Boundary/AMReX_LO_BCTYPES.H
    scal_matrix->setDomainBC(m_scalbc_lo, m_scalbc_hi);
    
    if(probtype==35){
        amrex::Print() << " warning tracer bcs are manually entered for probtype 35 in DiffusionOp.cpp " << std::endl;
        //fixme hard coded in
        Array<amrex::LinOpBCType,AMREX_SPACEDIM> tracbc_lo;
        Array<amrex::LinOpBCType,AMREX_SPACEDIM> tracbc_hi;
        tracbc_lo = {LinOpBCType::Periodic, LinOpBCType::Periodic, LinOpBCType::inhomogNeumann};
        tracbc_hi = {LinOpBCType::Periodic, LinOpBCType::Periodic, LinOpBCType::inhomogNeumann};
        
        // LinOpBCType Definitions are in amrex/Src/Boundary/AMReX_LO_BCTYPES.H
        scal_matrix->setDomainBC(tracbc_lo, tracbc_hi);
    }

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
// Solve the tensor equation for velocity
//
void DiffusionOp::diffuse_velocity(      Vector<std::unique_ptr<MultiFab>>& vel_in,
                                   const Vector<std::unique_ptr<MultiFab>>& ro_in,
                                   const Vector<std::unique_ptr<MultiFab>>& mu_in,
                                   Real dt)
{
    BL_PROFILE("DiffusionOp::diffuse_velocity");

    // Update the coefficients of the vel_matrix going into the solve based on the current state of the
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
    vel_matrix->setScalars(1.0, dt);

    int finest_level = amrcore->finestLevel();

    for(int lev = 0; lev <= finest_level; lev++)
    {
        // Compute the spatially varying b coefficients (on faces) to equal the apparent viscosity
        average_cellcenter_to_face(GetArrOfPtrs(b[lev]), *mu_in[lev], geom[lev]);
        for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
            b[lev][dir]->FillBoundary(geom[lev].periodicity());
        
        // This sets the coefficients
        vel_matrix->setACoeffs(lev, (*ro_in[lev]));
        vel_matrix->setShearViscosity  ( lev, GetArrOfConstPtrs(b[lev]));
#ifdef AMREX_USE_EB
        vel_matrix->setEBShearViscosity(lev, (*mu_in[lev]));
#endif
    }

    if(verbose > 0)
        amrex::Print() << "Diffusing velocity components all together..." << std::endl; 

    for(int lev = 0; lev <= finest_level; lev++)
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
        vel_matrix->setLevelBC(lev, GetVecOfConstPtrs(phi)[lev]);
    }

    MLMG solver(*vel_matrix);
    setSolverSettings(solver);

    // This ensures that ghost cells of sol are correctly filled when returned from the solver
    solver.setFinalFillBC(true);

    solver.solve(GetVecOfPtrs(phi), GetVecOfConstPtrs(rhs), mg_rtol, mg_atol);

    for(int lev = 0; lev <= finest_level; lev++)
    {
        phi[lev]->FillBoundary(geom[lev].periodicity());
        MultiFab::Copy(*vel_in[lev], *phi[lev], 0, 0, AMREX_SPACEDIM, 1);
    }
}

//
// Solve the equation for scalar diffusion
//
void DiffusionOp::diffuse_scalar(      Vector<std::unique_ptr<MultiFab>>& scal_in,
                                 const Vector<std::unique_ptr<MultiFab>>& ro_in,
                                 const Vector<std::unique_ptr<MultiFab>>& mu_in,
                                 Real dt)
{
    BL_PROFILE("DiffusionOp::diffuse_scalar");

    int finest_level = amrcore->finestLevel();

    // Update the coefficients of the scal_matrix going into the solve based on the current state of the
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
    scal_matrix->setScalars(1.0, dt);

    int ntrac = scal_in[0]->nComp();

    for(int lev = 0; lev <= finest_level; lev++)
    {
        average_cellcenter_to_face(GetArrOfPtrs(b[lev]), *mu_in[lev], geom[lev]);

        for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
            b[lev][dir]->FillBoundary(geom[lev].periodicity());
        
        // This sets the coefficients
        scal_matrix->setACoeffs (lev, (*ro_in[lev]));
        scal_matrix->setBCoeffs (lev, GetArrOfConstPtrs(b[lev]));
    }

    if(verbose > 0)
        amrex::Print() << "Diffusing tracers one at a time ..." << std::endl; 

    for(int lev = 0; lev <= finest_level; lev++)
    {
        // Zero these out just to have a clean start because they have 3 components 
        //      (due to re-use with velocity solve)
        phi[lev]->setVal(0.0);
        rhs[lev]->setVal(0.0);

        MultiFab::Copy((*rhs[lev]),(*scal_in[lev]), 0, 0, ntrac, 0);

        // Multiply rhs by rho  -- we are solving 
        //
        //      rho s_star = rho s_old + dt ( - div (rho u s) + rho div (mu (grad s)) )
        //
        for(int n = 0; n < ntrac; n++)
           MultiFab::Multiply((*rhs[lev]), (*ro_in[lev]), 0, n, 1, 0);

        MultiFab::Copy(*phi[lev],*scal_in[lev], 0, 0, ntrac, 1);
        scal_matrix->setLevelBC(lev, GetVecOfConstPtrs(phi)[lev]);
    }

    MLMG solver(*scal_matrix);
    setSolverSettings(solver);

    // This ensures that ghost cells of sol are correctly filled when returned from the solver
    solver.setFinalFillBC(true);

    solver.solve(GetVecOfPtrs(phi), GetVecOfConstPtrs(rhs), mg_rtol, mg_atol);

    for(int lev = 0; lev <= finest_level; lev++)
    {
        phi[lev]->FillBoundary(geom[lev].periodicity());
        MultiFab::Copy(*scal_in[lev], *phi[lev], 0, 0, ntrac, 1);
    }
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

    int finest_level = amrcore->finestLevel();

    Vector<std::unique_ptr<MultiFab> > divtau_aux(finest_level+1);
    for(int lev = 0; lev <= finest_level; lev++)
    {
#ifdef AMREX_USE_EB
       divtau_aux[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, nghost,
                                          MFInfo(), *(*ebfactory)[lev]));
#else
       divtau_aux[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, nghost, MFInfo()));
#endif
       divtau_aux[lev]->setVal(0.0);
    }
 
    // Whole domain
    Box domain(geom[0].Domain());
 
    // We want to return div (mu grad)) phi
    vel_matrix->setScalars(0.0, -1.0);
 
    // Compute the coefficients
    for (int lev = 0; lev <= finest_level; lev++)
    {
        average_cellcenter_to_face( GetArrOfPtrs(b[lev]), *mu_in[lev], geom[lev] );
 
        for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
             b[lev][dir]->FillBoundary(geom[lev].periodicity());
 
        vel_matrix->setShearViscosity  ( lev, GetArrOfConstPtrs(b[lev]));
#ifdef AMREX_USE_EB
        vel_matrix->setEBShearViscosity( lev, (*mu_in[lev]));
#endif
        vel_matrix->setLevelBC         ( lev, GetVecOfConstPtrs(vel_in)[lev] );
    }
 
    MLMG solver(*vel_matrix);
 
    solver.apply(GetVecOfPtrs(divtau_aux), GetVecOfPtrs(vel_in));
 
    for(int lev = 0; lev <= finest_level; lev++)
    {
#ifdef AMREX_USE_EB
       amrex::single_level_redistribute( lev, *divtau_aux[lev], *divtau_out[lev], 0, AMREX_SPACEDIM, geom);
#endif
 
       // Divide by density
       for (int n = 0; n < 3; n++)
           MultiFab::Divide( *divtau_out[lev], *ro_in[lev], 0, n, 1, 0 );
    }
}

void DiffusionOp::ComputeLapS(      Vector<std::unique_ptr<MultiFab>>& laps_out,
                              const Vector<std::unique_ptr<MultiFab>>& scal_in,
                              const Vector<std::unique_ptr<MultiFab>>& ro_in,
                              const Vector<Real> mu_s)
{
    BL_PROFILE("DiffusionOp::ComputeLapS");

    int finest_level = amrcore->finestLevel();

    int ntrac = scal_in[0]->nComp();

    Vector<std::unique_ptr<MultiFab> > laps_aux(finest_level+1);
    for(int lev = 0; lev <= finest_level; lev++)
    {
#ifdef AMREX_USE_EB
       laps_aux[lev].reset(new MultiFab(grids[lev], dmap[lev], ntrac, nghost,
                                          MFInfo(), *(*ebfactory)[lev]));
#else
       laps_aux[lev].reset(new MultiFab(grids[lev], dmap[lev], ntrac, nghost, MFInfo()));
#endif
       laps_aux[lev]->setVal(0.0);
    }
 
    // Whole domain
    Box domain(geom[0].Domain());
 
    // We want to return div (mu grad)) phi
    scal_matrix->setScalars(0.0, -1.0);
 
    // Compute the coefficients
    for (int lev = 0; lev <= finest_level; lev++)
    {
        for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
           for(int n = 0; n < ntrac; n++)
             b[lev][dir]->setVal(mu_s[n],n,1);
 
        scal_matrix->setBCoeffs  ( lev, GetArrOfConstPtrs(b[lev]));
        scal_matrix->setLevelBC  ( lev, GetVecOfConstPtrs(scal_in)[lev] );
    }
 
    MLMG solver(*scal_matrix);
 
    solver.apply(GetVecOfPtrs(laps_aux), GetVecOfPtrs(scal_in));
 
    for(int lev = 0; lev <= finest_level; lev++)
    {
#ifdef AMREX_USE_EB
       amrex::single_level_redistribute( lev, *laps_aux[lev], *laps_out[lev], 0, ntrac, geom);
#endif
 
       // Divide by density
       for(int n = 0; n < ntrac; n++)
          MultiFab::Divide( *laps_out[lev], *ro_in[lev], 0, n, 1, 0 );
    }
}

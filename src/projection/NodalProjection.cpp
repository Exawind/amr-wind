#include <NodalProjection.H>
#include <AMReX.H>
#include <incflo.H>
#include <AMReX_EBMultiFabUtil.H>

//
// Constructor:
// We set up everything which doesn't change between timesteps here
//
NodalProjection::NodalProjection(const incflo* a_incflo,
                                 Vector<std::unique_ptr<EBFArrayBoxFactory>>* _ebfactory,
                                 std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_bc_lo,
                                 std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_bc_hi)
{
    if (m_verbose > 0)
        amrex::Print() << "Constructing NodalProjection class" << std::endl;

    m_incflo = a_incflo;

    m_bc_lo = a_bc_lo;
    m_bc_hi = a_bc_hi;

    m_ok     = true;

    // Get inputs from ParmParse
    readParameters();

    // Actually do the setup work here
    setup(a_incflo, _ebfactory);
}

//
// Perform projection:
//
//     vel = vel - grad(phi)/ro
//
//  where phi is the solution of
//
//   div( sigma * grad(phi) ) = div(vel)
//
//  ro and vel are cell-centered
//
//  sigma = 1/ro is cell-centered as well
//
//  phi is node-centered
//
// If a_scale_factor is passed in, phi is return as phi/a_scale_factor
//
void
NodalProjection::project (      Vector< std::unique_ptr< amrex::MultiFab > >& a_vel,
                          const Vector< std::unique_ptr< amrex::MultiFab > >& a_ro,
                                Real a_time, Real a_scale_factor )
{
    AMREX_ALWAYS_ASSERT(m_ok);
    BL_PROFILE("NodalProjection::project");

    amrex::Print() << "Nodal Projection:" << std::endl;

    // Re-initialize phi before each projection
    for (int lev(0); lev < m_phi.size(); ++lev)
        m_phi[lev] -> setVal(0.0);

    // Compute RHS
    computeRHS(a_vel, a_time);

    // Print diagnostics
    amrex::Print() << " >> Before projection:" << std::endl;
    printInfo();

    // Compute and set matrix coefficients
    for (int lev(0); lev < a_ro.size(); ++lev)
    {
        // Compute the PPE coefficients = (1.0 / ro)
        m_sigma[lev] -> setVal(1.0);
        MultiFab::Divide(*m_sigma[lev],*a_ro[lev],0,0,1,0);

        // Set matrix coefficients
        m_matrix -> setSigma(lev, *m_sigma[lev]);
    }

    // Solve
    m_solver -> solve( GetVecOfPtrs(m_phi), GetVecOfConstPtrs(m_rhs), m_mg_rtol, m_mg_atol );

    // Get fluxes -- fluxes = - (1/ro)*grad(phi)
    m_solver -> getFluxes( GetVecOfPtrs(m_fluxes) );

    // Perform projection
    for (int lev(0); lev < m_phi.size(); ++lev)
    {
        // vel = vel + fluxes = vel - (scale_factor/rho) * grad(phi),
        MultiFab::Add( *a_vel[lev], *m_fluxes[lev], 0, 0, AMREX_SPACEDIM, 0);

        // Account for scale factor -- now fluxes = (1/rho) * grad(phi)
        m_fluxes[lev] -> mult(- 1.0/a_scale_factor, m_fluxes[lev]->nGrow() );

        // Finally we get rid of ro and MINUS so that m_fluxes = grad(phi)
        for (int n(0); n < AMREX_SPACEDIM; ++n)
            MultiFab::Multiply(*m_fluxes[lev], *a_ro[lev], 0, n, 1, m_fluxes[lev]->nGrow() );

        // Fill boundaries and apply scale factor to phi
        m_phi[lev] -> FillBoundary( geom[lev].periodicity());
        m_phi[lev] -> mult(1.0/a_scale_factor, m_fluxes[lev] -> nGrow());

    }

    // Compute RHS -- this is only needed to print out post projection values
    computeRHS(a_vel, a_time);

    // Print diagnostics
        amrex::Print() << " >> After projection:" << std::endl;
    printInfo();
}


//
// Read from input file
//
void
NodalProjection::readParameters ()
{
    ParmParse pp("projection");
    pp.query( "verbose"                , m_verbose );
    pp.query( "mg_verbose"             , m_mg_verbose );
    pp.query( "mg_cg_verbose"          , m_mg_cg_verbose );
    pp.query( "mg_maxiter"             , m_mg_maxiter );
    pp.query( "mg_cg_maxiter"          , m_mg_cg_maxiter );
    pp.query( "mg_rtol"                , m_mg_rtol );
    pp.query( "mg_atol"                , m_mg_atol );
    pp.query( "mg_max_coarsening_level", m_mg_max_coarsening_level );
    pp.query( "bottom_solver_type"     , m_bottom_solver_type );
}


//
// Setup object before solve
//
void
NodalProjection::setup(const incflo* a_incflo,
                       Vector<std::unique_ptr<EBFArrayBoxFactory>>* _ebfactory)
{
    BL_PROFILE("NodalProjection::setup");
    AMREX_ALWAYS_ASSERT(m_ok);

    // The incflo may change when we regrid so let's reset it here
    m_incflo = a_incflo;

    // The ebfactory changes when we regrid so we must pass it in here.
    ebfactory = _ebfactory;

    geom  = m_incflo->Geom();
    grids = m_incflo->boxArray();
    dmap  = m_incflo->DistributionMap();

    // Set number of levels
    int nlev( grids.size() );

    // Resize member data if necessary
    if ( nlev != m_phi.size() )
    {
        m_phi.resize(nlev);
        m_fluxes.resize(nlev);
        m_sigma.resize(nlev);
        m_rhs.resize(nlev);
    }

    // Regrid if necessary
    int nghost(1);      // We use 1 ghost node only -- it should be enough

    bool need_regrid(false);  // if BA and DM changed on any level, we need
                              // to update the matrix and the solver as well

    for (int lev(0); lev < nlev; ++lev )
    {
        const auto& eb = *(m_incflo -> ebfactory[lev]);

        if ( (m_phi[lev] == nullptr)                 ||
             (m_phi[lev] -> boxArray()        != grids[lev]) ||
             (m_phi[lev] -> DistributionMap() != dmap[lev])  )
        {
            // Cell-centered data
            m_fluxes[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, nghost, MFInfo(), eb));
             m_sigma[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost, MFInfo(), eb));

            // Node-centered data
            const auto& ba_nd = amrex::convert(grids[lev], IntVect{1,1,1});
            m_phi[lev].reset(new MultiFab(ba_nd, dmap[lev], 1, nghost, MFInfo(), eb));
            m_rhs[lev].reset(new MultiFab(ba_nd, dmap[lev], 1, nghost, MFInfo(), eb));

            need_regrid = true;
        }
    }

    // Setup matrix and solver
    if ( (m_matrix == nullptr) || need_regrid )
    {
        //
        // Setup Matrix
        //
        LPInfo                       info;
        info.setMaxCoarseningLevel(m_mg_max_coarsening_level);
        m_matrix.reset(new MLNodeLaplacian(geom, grids, dmap, info,
                                           GetVecOfConstPtrs(m_incflo->ebfactory)));

        m_matrix->setGaussSeidel(true);
        m_matrix->setHarmonicAverage(false);
        m_matrix->setDomainBC(m_bc_lo, m_bc_hi);

        //
        // Setup solver
        //
        m_solver.reset(new MLMG(*m_matrix));

        m_solver->setMaxIter(m_mg_maxiter);
        m_solver->setVerbose(m_mg_verbose);
        m_solver->setCGVerbose(m_mg_cg_verbose);
        m_solver->setCGMaxIter(m_mg_cg_maxiter);

        if (m_bottom_solver_type == "smoother")
        {
            m_solver->setBottomSolver(MLMG::BottomSolver::smoother);
        }
        else if (m_bottom_solver_type == "bicg")
        {
            m_solver->setBottomSolver(MLMG::BottomSolver::bicgstab);
        }
        else if (m_bottom_solver_type == "cg")
        {
            m_solver->setBottomSolver(MLMG::BottomSolver::cg);
        }
        else if (m_bottom_solver_type == "bicgcg")
        {
            m_solver->setBottomSolver(MLMG::BottomSolver::bicgcg);
        }
        else if (m_bottom_solver_type == "cgbicg")
        {
            m_solver->setBottomSolver(MLMG::BottomSolver::cgbicg);
        }
        else if (m_bottom_solver_type == "hypre")
        {
#ifdef AMREX_USE_HYPRE
            m_solver->setBottomSolver(MLMG::BottomSolver::hypre);
#else
            amrex::Abort("AMReX was not built with HYPRE support");
#endif
        }
    }
}

//
// Return DivU for diagnostics
//
void
NodalProjection::getDivU (Vector< std::unique_ptr< amrex::MultiFab > >& divu,
                          Vector< std::unique_ptr< amrex::MultiFab > >& a_vel,
                          Real a_time ) 
{
    AMREX_ALWAYS_ASSERT(m_ok);
    m_matrix -> compDivergence( GetVecOfPtrs(divu),  GetVecOfPtrs(a_vel));
}

//
// Compute RHS: div(u) (later this may have a specified S as in div(u) = S)
//

void
NodalProjection::computeRHS ( Vector< std::unique_ptr< amrex::MultiFab > >& a_vel,
                              Real a_time )
{
    AMREX_ALWAYS_ASSERT(m_ok);
    BL_PROFILE("NodalProjection::computeRHS");

    int extrap_dir_bcs(0);
    m_incflo -> incflo_set_velocity_bcs(a_time, a_vel, extrap_dir_bcs);

    // Compute div(eu)
    m_matrix -> compRHS( GetVecOfPtrs(m_rhs),  GetVecOfPtrs(a_vel), {}, {} );
}


void
NodalProjection::printInfo ()
{
    for (int lev(0); lev < m_rhs.size(); ++lev)
    {
        amrex::Print() << "  * On lev " << lev
                       << " max(abs(divu)) = "
                       << m_rhs[lev]->norm0(0,0,false,true)
                       << std::endl;
    }
}

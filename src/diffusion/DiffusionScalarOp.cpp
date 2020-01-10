#include <incflo.H>
#include <DiffusionScalarOp.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

DiffusionScalarOp::DiffusionScalarOp (incflo* a_incflo)
    : m_incflo(a_incflo)
{
    readParameters();

    int finest_level = m_incflo->finestLevel();

    LPInfo info;
    info.setMaxCoarseningLevel(m_mg_max_coarsening_level);
#ifdef AMREX_USE_EB
    if (!m_incflo->EBFactory(0).isAllRegular())
    {
        Vector<EBFArrayBoxFactory const*> ebfact;
        for (int lev = 0; lev <= finest_level; ++lev) {
            ebfact.push_back(&(m_incflo->EBFactory(lev)));
        }
        m_eb_op.reset(new MLEBABecLap(m_incflo->Geom(0,finest_level),
                                      m_incflo->boxArray(0,finest_level),
                                      m_incflo->DistributionMap(0,finest_level),
                                      info, ebfact));
        m_eb_op->setMaxOrder(m_mg_maxorder);
        m_eb_op->setDomainBC(m_incflo->get_diffuse_scalar_bc(Orientation::low),
                             m_incflo->get_diffuse_scalar_bc(Orientation::high));
    }
    else
#endif
    {
        m_reg_op.reset(new MLABecLaplacian(m_incflo->Geom(0,m_incflo->finestLevel()),
                                           m_incflo->boxArray(0,m_incflo->finestLevel()),
                                           m_incflo->DistributionMap(0,m_incflo->finestLevel()),
                                           info));
        m_reg_op->setMaxOrder(m_mg_maxorder);
        m_reg_op->setDomainBC(m_incflo->get_diffuse_scalar_bc(Orientation::low),
                              m_incflo->get_diffuse_scalar_bc(Orientation::high));
    }
}

void
DiffusionScalarOp::readParameters ()
{
    ParmParse pp("diffusion");

    pp.query("verbose", m_verbose);
    pp.query("mg_verbose", m_mg_verbose);
    pp.query("mg_cg_verbose", m_mg_cg_verbose);
    pp.query("mg_max_iter", m_mg_max_iter);
    pp.query("mg_cg_maxiter", m_mg_cg_maxiter);
    pp.query("mg_max_fmg_iter", m_mg_max_fmg_iter);
    pp.query("mg_max_coarsening_level", m_mg_max_coarsening_level);
    pp.query("mg_maxorder", m_mg_maxorder);
    pp.query("mg_rtol", m_mg_rtol);
    pp.query("mg_atol", m_mg_atol);
    pp.query("bottom_solver_type", m_bottom_solver_type);
}

void
DiffusionScalarOp::diffuse_scalar (Vector<MultiFab*> const& tracer,
                                   Vector<MultiFab*> const& density,
                                   Real t, Real dt)
{
    //
    //      alpha a - beta div ( b grad )   <--->   rho - dt div ( mu grad )
    //
    // So the constants and variable coefficients are:
    //
    //      alpha: 1
    //      beta: dt
    //      a: rho
    //      b: mu

    if (m_verbose > 0) {
        amrex::Print() << "Diffusing scalars one at a time ..." << std::endl;
    }

    const int finest_level = m_incflo->finestLevel();

    Vector<MultiFab> rhs(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        rhs[lev].define(tracer[lev]->boxArray(), tracer[lev]->DistributionMap(), 1, 0);
    }

#ifdef AMREX_USE_EB
    if (m_eb_op)
    {
        m_eb_op->setScalars(1.0, dt);
        for (int lev = 0; lev <= finest_level; ++lev) {
            m_eb_op->setACoeffs(lev, *density[lev]);
        }
    }
    else
#endif
    {
        m_reg_op->setScalars(1.0, dt);
        for (int lev = 0; lev <= finest_level; ++lev) {
            m_reg_op->setACoeffs(lev, *density[lev]);
        }
    }

    for (int comp = 0; comp < tracer[0]->nComp(); ++comp)
    {
#ifdef AMREX_USE_EB
        if (m_eb_op)
        {
            for (int lev = 0; lev <= finest_level; ++lev) {
                m_eb_op->setBCoeffs(lev, m_incflo->mu_s[comp]);
            }
        }
        else
#endif
        {
            for (int lev = 0; lev <= finest_level; ++lev) {
                m_reg_op->setBCoeffs(lev, m_incflo->mu_s[comp]);
            }
        }

        Vector<MultiFab> phi;
        for (int lev = 0; lev <= finest_level; ++lev) {
            phi.emplace_back(*tracer[lev], amrex::make_alias, comp, 1);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(rhs[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                Box const& bx = mfi.tilebox();
                Array4<Real> const& rhs_a = rhs[lev].array(mfi);
                Array4<Real const> const& tra_a = tracer[lev]->const_array(mfi,comp);
                Array4<Real const> const& rho_a = density[lev]->const_array(mfi);
                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    rhs_a(i,j,k) = rho_a(i,j,k) * tra_a(i,j,k);
                });
            }

#ifdef AMREX_USE_EB
            if (m_eb_op) {
                m_eb_op->setLevelBC(lev, tracer[lev]);
            } else
#endif
            {
                m_reg_op->setLevelBC(lev, tracer[lev]);
            }
        }

#ifdef AMREX_USE_EB
        MLMG mlmg(m_eb_op ? static_cast<MLLinOp&>(*m_eb_op) : static_cast<MLLinOp&>(*m_reg_op));
#else
        MLMG mlmg(*m_reg_op);
#endif

        // The default bottom solver is BiCG
        if (m_bottom_solver_type == "smoother")
        {
            mlmg.setBottomSolver(MLMG::BottomSolver::smoother);
        }
        else if (m_bottom_solver_type == "hypre")
        {
            mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
        }
        // Maximum iterations for MultiGrid / ConjugateGradients
        mlmg.setMaxIter(m_mg_max_iter);
        mlmg.setMaxFmgIter(m_mg_max_fmg_iter);
        mlmg.setCGMaxIter(m_mg_cg_maxiter);
        
        // Verbosity for MultiGrid / ConjugateGradients
        mlmg.setVerbose(m_mg_verbose);
        mlmg.setCGVerbose(m_mg_cg_verbose);

        amrex::EB_set_covered(*tracer[0], 0.0);
        VisMF::Write(*tracer[0], "tra0");

        amrex::system::verbose = 10;
        amrex::OutStream().precision(17);

        mlmg.solve(GetVecOfPtrs(phi), GetVecOfConstPtrs(rhs), m_mg_rtol, m_mg_atol);

        amrex::EB_set_covered(*tracer[0], 0.0);
        VisMF::Write(*tracer[0], "tra1");

        amrex::Abort("xxxxxx after scalar solve");
    }
}

#include <incflo.H>
#include <DiffusionScalarOp.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

DiffusionScalarOp::DiffusionScalarOp (incflo* a_incflo)
    : m_incflo(a_incflo)
{
    readParameters();

    LPInfo info_solve;
    info_solve.setMaxCoarseningLevel(m_mg_max_coarsening_level);
    LPInfo info_apply;
    info_apply.setMaxCoarseningLevel(0);

    m_reg_solve_op.reset(new MLABecLaplacian(m_incflo->Geom(0,m_incflo->finestLevel()),
                                             m_incflo->boxArray(0,m_incflo->finestLevel()),
                                             m_incflo->DistributionMap(0,m_incflo->finestLevel()),
                                             info_solve));
    m_reg_solve_op->setMaxOrder(m_mg_maxorder);
    m_reg_solve_op->setDomainBC(m_incflo->get_diffuse_scalar_bc(Orientation::low),
                                m_incflo->get_diffuse_scalar_bc(Orientation::high));
    if (m_incflo->need_divtau()) {
        m_reg_apply_op.reset(new MLABecLaplacian(m_incflo->Geom(0,m_incflo->finestLevel()),
                                                 m_incflo->boxArray(0,m_incflo->finestLevel()),
                                                 m_incflo->DistributionMap(0,m_incflo->finestLevel()),
                                                 info_apply));
        m_reg_apply_op->setMaxOrder(m_mg_maxorder);
        m_reg_apply_op->setDomainBC(m_incflo->get_diffuse_scalar_bc(Orientation::low),
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
                                   Vector<MultiFab const*> const& eta,
                                   Real dt)
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
    
    m_reg_solve_op->setScalars(1.0, dt);
    for (int lev = 0; lev <= finest_level; ++lev) {
        m_reg_solve_op->setACoeffs(lev, *density[lev]);
    }
    

    for (int comp = 0; comp < tracer[0]->nComp(); ++comp)
    {

        Vector<MultiFab> phi;
        for (int lev = 0; lev <= finest_level; ++lev) {
            phi.emplace_back(*tracer[lev], amrex::make_alias, comp, 1);

            Array<MultiFab,AMREX_SPACEDIM> b = m_incflo->average_tracer_eta_to_faces(lev, comp, *eta[lev]);
            m_reg_solve_op->setBCoeffs(lev, GetArrOfConstPtrs(b));
            
            MultiFab bc(phi[lev].boxArray(),phi[lev].DistributionMap(),1, 1, MFInfo(),phi[lev].Factory());

            // copy tracer into bc in case there is a periodic bc
            // fixme this seems inefficient
            MultiFab::Copy(bc, phi[lev], 0, 0, 1, 1);
            m_incflo->heat_flux_model_bc(lev, comp, bc);
            m_reg_solve_op->setLevelBC(lev, &bc);
            
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
        }

        MLMG mlmg(*m_reg_solve_op);

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

        mlmg.solve(GetVecOfPtrs(phi), GetVecOfConstPtrs(rhs), m_mg_rtol, m_mg_atol);
    }
}

void DiffusionScalarOp::compute_laps (Vector<MultiFab*> const& a_laps,
                                      Vector<MultiFab const*> const& a_tracer,
                                      Vector<MultiFab const*> const& a_density,
                                      Vector<MultiFab const*> const& a_eta)
{
    BL_PROFILE("DiffusionScalarOp::compute_laps");

    int finest_level = m_incflo->finestLevel();

    Vector<MultiFab> tracer(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        tracer[lev].define(a_tracer[lev]->boxArray(),
                           a_tracer[lev]->DistributionMap(),
                           m_incflo->m_ntrac, 1, MFInfo(),
                           a_tracer[lev]->Factory());
        MultiFab::Copy(tracer[lev], *a_tracer[lev], 0, 0, m_incflo->m_ntrac, 1);
    }
    
    // We want to return div (mu grad)) phi
    m_reg_apply_op->setScalars(0.0, -1.0);

    // This should have no effect since the first scalar is 0
    for (int lev = 0; lev <= finest_level; ++lev) {
        m_reg_apply_op->setACoeffs(lev, *a_density[lev]);
    }

    for (int comp = 0; comp < m_incflo->m_ntrac; ++comp) {
        Vector<MultiFab> laps_comp;
        Vector<MultiFab> tracer_comp;
        for (int lev = 0; lev <= finest_level; ++lev) {
            laps_comp.emplace_back(*a_laps[lev],amrex::make_alias,comp,1);
            tracer_comp.emplace_back(tracer[lev],amrex::make_alias,comp,1);
            Array<MultiFab,AMREX_SPACEDIM> b = m_incflo->average_tracer_eta_to_faces(lev, comp, *a_eta[lev]);
            m_reg_apply_op->setBCoeffs(lev, GetArrOfConstPtrs(b));
            
            MultiFab bc(tracer_comp[lev].boxArray(), tracer_comp[lev].DistributionMap(), 1, 1, MFInfo(), tracer_comp[lev].Factory());
            // copy tracer into bc in case there is a periodic bc
            // fixme this seems inefficient
            MultiFab::Copy(bc, tracer_comp[lev], 0, 0, 1, 1);

            m_incflo->heat_flux_model_bc(lev, comp, bc);
            m_reg_apply_op->setLevelBC(lev, &bc);
        }

        MLMG mlmg(*m_reg_apply_op);
        mlmg.apply(GetVecOfPtrs(laps_comp), GetVecOfPtrs(tracer_comp));
    }
    
}

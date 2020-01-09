#include <incflo.H>
#include <DiffusionTensorOp.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

DiffusionTensorOp::DiffusionTensorOp (incflo* a_incflo)
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
        m_eb_op.reset(new MLEBTensorOp(m_incflo->Geom(0,finest_level),
                                       m_incflo->boxArray(0,finest_level),
                                       m_incflo->DistributionMap(0,finest_level),
                                       info, ebfact));
        m_eb_op->setMaxOrder(m_mg_maxorder);
        m_eb_op->setDomainBC(m_incflo->get_diffuse_tensor_bc(Orientation::low),
                             m_incflo->get_diffuse_tensor_bc(Orientation::high));
    }
    else
#endif
    {
        m_reg_op.reset(new MLTensorOp(m_incflo->Geom(0,m_incflo->finestLevel()),
                                      m_incflo->boxArray(0,m_incflo->finestLevel()),
                                      m_incflo->DistributionMap(0,m_incflo->finestLevel()),
                                      info));
        m_reg_op->setMaxOrder(m_mg_maxorder);
        m_reg_op->setDomainBC(m_incflo->get_diffuse_tensor_bc(Orientation::low),
                              m_incflo->get_diffuse_tensor_bc(Orientation::high));
    }
}

void
DiffusionTensorOp::readParameters ()
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
DiffusionTensorOp::diffuse_velocity (Vector<MultiFab*> const& velocity,
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
        amrex::Print() << "Diffusing velocity components all together..." << std::endl;
    }

    const int finest_level = m_incflo->finestLevel();

#ifdef AMREX_USE_EB
    if (m_eb_op)
    {
        m_eb_op->setScalars(1.0, dt);
        for (int lev = 0; lev <= finest_level; ++lev) {
            m_eb_op->setACoeffs(lev, *density[lev]);
            m_eb_op->setShearViscosity(lev, m_incflo->mu);
            m_eb_op->setEBShearViscosity(lev, m_incflo->mu);
        }
    }
    else
#endif
    {
        m_reg_op->setScalars(1.0, dt);
        for (int lev = 0; lev <= finest_level; ++lev) {
            m_reg_op->setACoeffs(lev, *density[lev]);
            m_reg_op->setShearViscosity(lev, m_incflo->mu);
        }
    }

    amrex::Abort("xxxxx so far so good in diffuse_velocity");
}

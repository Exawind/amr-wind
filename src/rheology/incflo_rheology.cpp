#include <incflo.H>
#include <derive_K.H>

using namespace amrex;

namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real expterm (amrex::Real nu) noexcept
{
    return (nu < 1.e-9) ? (1.0-0.5*nu+nu*nu*(1.0/6.0)-(nu*nu*nu)*(1./24.))
                        : -std::expm1(-nu)/nu;
}

struct NonNewtonianViscosity
{
    incflo::FluidModel fluid_model;
    amrex::Real mu, n_flow, tau_0, eta_0, papa_reg;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real operator() (amrex::Real sr) const noexcept {
        switch (fluid_model)
        {
        case incflo::FluidModel::powerlaw:
        {
            return mu * std::pow(sr,n_flow-1.0);
        }
        case incflo::FluidModel::Bingham:
        {
            return mu + tau_0 * expterm(sr/papa_reg) / papa_reg;
        }
        case incflo::FluidModel::HerschelBulkley:
        {
            return (mu*std::pow(sr,n_flow)+tau_0)*expterm(sr/papa_reg)/papa_reg;
        }
        case incflo::FluidModel::deSouzaMendesDutra:
        {
            return (mu*std::pow(sr,n_flow)+tau_0)*expterm(sr*(eta_0/tau_0))*(eta_0/tau_0);
        }
        default:
        {
            return mu;
        }
        };
    }
};

}

void incflo::compute_viscosity (Vector<MultiFab*> const& vel_eta,
                                Vector<MultiFab*> const& tra_eta,
                                Vector<MultiFab const*> const& rho,
                                Vector<MultiFab const*> const& vel,
                                Vector<MultiFab const*> const& tra,
                                Real time, int nghost)
{
    if (m_fluid_model == FluidModel::Newtonian)
    {
        for (auto mf : vel_eta) {
            mf->setVal(m_mu, 0, 1, nghost);
        }
    }
    else
    {
        NonNewtonianViscosity non_newtonian_viscosity;
        non_newtonian_viscosity.fluid_model = m_fluid_model;
        non_newtonian_viscosity.mu = m_mu;
        non_newtonian_viscosity.n_flow = m_n_0;
        non_newtonian_viscosity.tau_0 = m_tau_0;
        non_newtonian_viscosity.eta_0 = m_eta_0;
        non_newtonian_viscosity.papa_reg = m_papa_reg;

        for (int lev = 0; lev <= finest_level; ++lev) {
#ifdef AMREX_USE_EB
            auto const& fact = EBFactory(lev);
            auto const& flags = fact.getMultiEBCellFlagFab();
#endif

            Real idx = 1.0 / geom[lev].CellSize(0);
            Real idy = 1.0 / geom[lev].CellSize(1);
            Real idz = 1.0 / geom[lev].CellSize(2);

#ifdef _OPENMP
#pragma omp parallel omp if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*vel_eta[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box const& bx = mfi.growntilebox(nghost);
                Array4<Real> const& eta_arr = vel_eta[lev]->array(mfi);
                Array4<Real const> const& vel_arr = vel[lev]->const_array(mfi);
#ifdef AMREX_USE_EB
                auto const& flag_fab = flags[mfi];
                auto typ = flag_fab.getType(bx);
                if (typ == FabType::covered)
                {
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        eta_arr(i,j,k) = 0.0;
                    });
                }
                else if (typ == FabType::singlevalued)
                {
                    auto const& flag_arr = flag_fab.const_array();
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        Real sr = incflo_strainrate_eb(i,j,k,idx,idy,idz,vel_arr,flag_arr(i,j,k));
                        eta_arr(i,j,k) = non_newtonian_viscosity(sr);
                    });
                }
                else
#endif
                {
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        Real sr = incflo_strainrate(i,j,k,idx,idy,idz,vel_arr);
                        eta_arr(i,j,k) = non_newtonian_viscosity(sr);
                    });
                }
            }
        }
    }

    for (auto mf : tra_eta) {
        for (int n = 0; n < m_ntrac; ++n) {
            mf->setVal(m_mu_s[n], n, 1, nghost);
        }
    }
}

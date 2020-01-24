#include <incflo.H>

using namespace amrex;

void incflo::compute_viscosity (Vector<MultiFab*> const& vel_eta,
                                Vector<MultiFab*> const& tra_eta,
                                Vector<MultiFab const*> const& rho,
                                Vector<MultiFab const*> const& vel,
                                Vector<MultiFab const*> const& tra,
                                Real time, int nghost)
{
    if (m_fluid_model == "newtonian")
    {
        for (auto mf : vel_eta) {
            mf->setVal(m_mu, 0, 1, nghost);
        }
    }
    else
    {
        amrex::Abort("xxxxx TODO");
    }

    for (auto mf : tra_eta) {
        for (int n = 0; n < m_ntrac; ++n) {
            mf->setVal(m_mu_s[n], n, 1, nghost);
        }
    }
}

#include "amr-wind/mesh_mapping_models/ConstantScaling.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

ConstantScaling::ConstantScaling(const CFDSim& sim)
    : m_mesh_scale_fac(sim.repo().get_field("mesh_scaling_factor"))
{
    amrex::ParmParse pp("ConstantScaling");
    pp.getarr("scaling_factor",m_fac,0,AMREX_SPACEDIM);
}

/** Construct the mesh mapping field
 */
void ConstantScaling::create_map(int finest_level)
{
    amrex::Real fac_x = m_fac[0];
    amrex::Real fac_y = m_fac[1];
    amrex::Real fac_z = m_fac[2];

    for (int lev = 0; lev <= finest_level; ++lev) {
        for (amrex::MFIter mfi(m_mesh_scale_fac(lev)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.validbox();
            amrex::Array4<amrex::Real> const& scale_fac = m_mesh_scale_fac(lev).array(mfi);

            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    scale_fac(i, j, k, 0) = fac_x;
                    scale_fac(i, j, k, 1) = fac_y;
                    scale_fac(i, j, k, 2) = fac_z;
                });
        }
    }
}

} // namespace amr_wind

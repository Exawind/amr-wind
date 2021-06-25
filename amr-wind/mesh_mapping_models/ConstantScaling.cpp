#include "amr-wind/mesh_mapping_models/ConstantScaling.H"
#include "amr-wind/CFDSim.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {

ConstantScaling::ConstantScaling(const CFDSim& sim)
    : m_mesh_scale_fac_cc(sim.repo().get_field("mesh_scaling_factor_cc"))
    , m_mesh_scale_fac_nd(sim.repo().get_field("mesh_scaling_factor_nd"))
{
    amrex::ParmParse pp("ConstantScaling");
    if(pp.contains("ConstantScaling")) {
        pp.getarr("scaling_factor",m_fac,0,AMREX_SPACEDIM);
    }
}

/** Construct the mesh mapping field
 */
void ConstantScaling::create_map(
    int lev, const amrex::Geometry& geom)
{
    amrex::Real fac_x = m_fac[0];
    amrex::Real fac_y = m_fac[1];
    amrex::Real fac_z = m_fac[2];

    for (amrex::MFIter mfi(m_mesh_scale_fac_cc(lev)); mfi.isValid(); ++mfi) {

        const auto& bx = mfi.growntilebox();
        amrex::Array4<amrex::Real> const& scale_fac_cc = m_mesh_scale_fac_cc(lev).array(mfi);
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                scale_fac_cc(i, j, k, 0) = fac_x;
                scale_fac_cc(i, j, k, 1) = fac_y;
                scale_fac_cc(i, j, k, 2) = fac_z;
            });

        const auto& nbx = mfi.grownnodaltilebox();
        amrex::Array4<amrex::Real> const& scale_fac_nd = m_mesh_scale_fac_nd(lev).array(mfi);
        amrex::ParallelFor(
            nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                scale_fac_nd(i, j, k, 0) = fac_x;
                scale_fac_nd(i, j, k, 1) = fac_y;
                scale_fac_nd(i, j, k, 2) = fac_z;
            });
    }
}

} // namespace amr_wind

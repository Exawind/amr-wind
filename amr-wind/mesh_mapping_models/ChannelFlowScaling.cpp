#include "amr-wind/mesh_mapping_models/ChannelFlowScaling.H"
#include "amr-wind/CFDSim.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {

ChannelFlowScaling::ChannelFlowScaling(const CFDSim& sim)
    : m_mesh_scale_fac_cc(sim.repo().get_field("mesh_scaling_factor_cc"))
    , m_mesh_scale_fac_nd(sim.repo().get_field("mesh_scaling_factor_nd"))
{
    amrex::ParmParse pp("ChannelFlowScaling");
    pp.queryarr("alpha",m_alpha,0,AMREX_SPACEDIM);
    pp.queryarr("beta",m_beta,0,AMREX_SPACEDIM);
    pp.queryarr("do_map",m_map,0,AMREX_SPACEDIM);
}

/** Construct the mesh mapping field
 */
void ChannelFlowScaling::create_map(
    int lev, const amrex::Geometry& geom)
{
    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();
    const auto& prob_hi = geom.ProbLoArray();
    amrex::Vector<amrex::Real> ht{{
        prob_hi[0] - prob_lo[0],
        prob_hi[1] - prob_lo[1],
        prob_hi[2] - prob_lo[2]
    }};

    for (amrex::MFIter mfi(m_mesh_scale_fac_cc(lev)); mfi.isValid(); ++mfi) {

        const auto& bx = mfi.growntilebox();
        amrex::Array4<amrex::Real> const& scale_fac_cc = m_mesh_scale_fac_cc(lev).array(mfi);
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real x = prob_lo[0] + (i+0.5) * dx[0];
                amrex::Real y = prob_lo[1] + (j+0.5) * dx[1];
                amrex::Real z = prob_lo[2] + (k+0.5) * dx[2];

                amrex::Real num_x = ht[0] * (m_beta[0] + 2*m_alpha[0])
                                  * std::pow(((m_beta[0]+1)/(m_beta[0]-1)),((x-m_alpha[0])/(1-m_alpha[0])))
                                  - m_beta[0] + 2*m_alpha[0];
                amrex::Real num_y = ht[1] * (m_beta[1] + 2*m_alpha[1])
                                  * std::pow(((m_beta[1]+1)/(m_beta[1]-1)),((y-m_alpha[1])/(1-m_alpha[1])))
                                  - m_beta[1] + 2*m_alpha[1];
                amrex::Real num_z = ht[2] * (m_beta[2] + 2*m_alpha[2])
                                  * std::pow(((m_beta[2]+1)/(m_beta[2]-1)),((z-m_alpha[2])/(1-m_alpha[2])))
                                  - m_beta[2] + 2*m_alpha[2];

                amrex::Real den_x = (2*m_alpha[0] + 1) * (1
                                  + std::pow(((m_beta[0]+1)/(m_beta[0]-1)),((x-m_alpha[0])/(1-m_alpha[0]))));
                amrex::Real den_y = (2*m_alpha[1] + 1) * (1
                                  + std::pow(((m_beta[1]+1)/(m_beta[1]-1)),((y-m_alpha[1])/(1-m_alpha[1]))));
                amrex::Real den_z = (2*m_alpha[2] + 1) * (1
                                  + std::pow(((m_beta[2]+1)/(m_beta[2]-1)),((z-m_alpha[2])/(1-m_alpha[2]))));

                scale_fac_cc(i, j, k, 0) = m_map[0] ? num_x/den_x/dx[0] : 1.0;
                scale_fac_cc(i, j, k, 1) = m_map[1] ? num_y/den_y/dx[1] : 1.0;
                scale_fac_cc(i, j, k, 2) = m_map[2] ? num_z/den_z/dx[2] : 1.0;
            });

        const auto& nbx = mfi.grownnodaltilebox();
        amrex::Array4<amrex::Real> const& scale_fac_nd = m_mesh_scale_fac_nd(lev).array(mfi);
        amrex::ParallelFor(
            nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real x = prob_lo[0] + i * dx[0];
                amrex::Real y = prob_lo[1] + j * dx[1];
                amrex::Real z = prob_lo[2] + k * dx[2];

                amrex::Real num_x = ht[0] * (m_beta[0] + 2*m_alpha[0])
                                  * std::pow(((m_beta[0]+1)/(m_beta[0]-1)),((x-m_alpha[0])/(1-m_alpha[0])))
                                  - m_beta[0] + 2*m_alpha[0];
                amrex::Real num_y = ht[1] * (m_beta[1] + 2*m_alpha[1])
                                  * std::pow(((m_beta[1]+1)/(m_beta[1]-1)),((y-m_alpha[1])/(1-m_alpha[1])))
                                  - m_beta[1] + 2*m_alpha[1];
                amrex::Real num_z = ht[2] * (m_beta[2] + 2*m_alpha[2])
                                  * std::pow(((m_beta[2]+1)/(m_beta[2]-1)),((z-m_alpha[2])/(1-m_alpha[2])))
                                  - m_beta[2] + 2*m_alpha[2];

                amrex::Real den_x = (2*m_alpha[0] + 1) * (1
                                  + std::pow(((m_beta[0]+1)/(m_beta[0]-1)),((x-m_alpha[0])/(1-m_alpha[0]))));
                amrex::Real den_y = (2*m_alpha[1] + 1) * (1
                                  + std::pow(((m_beta[1]+1)/(m_beta[1]-1)),((y-m_alpha[1])/(1-m_alpha[1]))));
                amrex::Real den_z = (2*m_alpha[2] + 1) * (1
                                  + std::pow(((m_beta[2]+1)/(m_beta[2]-1)),((z-m_alpha[2])/(1-m_alpha[2]))));

                scale_fac_nd(i, j, k, 0) = m_map[0] ? num_x/den_x/dx[0] : 1.0;
                scale_fac_nd(i, j, k, 1) = m_map[1] ? num_y/den_y/dx[1] : 1.0;
                scale_fac_nd(i, j, k, 2) = m_map[2] ? num_z/den_z/dx[2] : 1.0;
            });
    }
}

} // namespace amr_wind

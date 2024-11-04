
#include "amr-wind/equation_systems/temperature/source_terms/TempSpongeForcing.H"
#include "amr-wind/utilities/IOManager.H"
#include "amr-wind/utilities/linear_interpolation.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include "AMReX_Random.H"

namespace amr_wind::pde::temperature {

TempSpongeForcing::TempSpongeForcing(const CFDSim& sim)
    : m_mesh(sim.mesh()), m_temperature(sim.repo().get_field("temperature"))
{
    amrex::ParmParse pp_abl("ABL");
    //! Temperature variation as a function of height
    pp_abl.query("meso_sponge_start", m_sponge_start);
    pp_abl.getarr("temperature_heights", m_theta_heights);
    pp_abl.getarr("temperature_values", m_theta_values);
    AMREX_ALWAYS_ASSERT(m_theta_heights.size() == m_theta_values.size());
    const int num_theta_values = static_cast<int>(m_theta_heights.size());
    m_theta_heights_d.resize(num_theta_values);
    m_theta_values_d.resize(num_theta_values);
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_theta_heights.begin(),
        m_theta_heights.end(), m_theta_heights_d.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_theta_values.begin(), m_theta_values.end(),
        m_theta_values_d.begin());
}

TempSpongeForcing::~TempSpongeForcing() = default;

void TempSpongeForcing::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState fstate,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto& geom = m_mesh.Geom(lev);
    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();
    const auto& prob_hi = geom.ProbHiArray();
    const auto& temperature =
        m_temperature.state(field_impl::dof_state(fstate))(lev).const_array(
            mfi);
    const amrex::Real sponge_start = m_sponge_start;
    const auto vsize = m_theta_heights_d.size();
    const auto* theta_heights_d = m_theta_heights_d.data();
    const auto* theta_values_d = m_theta_values_d.data();
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
        const amrex::Real zi =
            std::max((z - sponge_start) / (prob_hi[2] - sponge_start), 0.0);
        amrex::Real ref_temp = temperature(i, j, k);
        if (zi > 0) {
            const auto idx = interp::bisection_search(
                theta_heights_d, theta_heights_d + vsize, z);
            ref_temp = (vsize > 0)
                           ? interp::linear_impl(
                                 theta_heights_d, theta_values_d, z, idx)
                           : temperature(i, j, k);
        }
        src_term(i, j, k, 0) -= zi * zi * (temperature(i, j, k) - ref_temp);
    });
}

} // namespace amr_wind::pde::temperature

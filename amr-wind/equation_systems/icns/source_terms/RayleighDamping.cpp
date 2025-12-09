#include "amr-wind/equation_systems/icns/source_terms/RayleighDamping.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/trig_ops.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"

namespace amr_wind::pde::icns {

RayleighDamping::RayleighDamping(const CFDSim& sim)
    : m_mesh(sim.mesh()), m_velocity(sim.repo().get_field("velocity"))
{
    // Read the Rayleigh Damping Layer parameters
    amrex::ParmParse pp("RayleighDamping");
    pp.get("time_scale", m_tau);
    // Length where damping coefficient depends on spatial position
    // In sloped region, coefficient goes from 1 to 0
    pp.get("length_sloped_damping", m_dRD);
    // Length where damping coefficient is set to 1
    pp.get("length_complete_damping", m_dFull);
    // Total damping length is m_dRD + m_dFull. Total length is not read in.
    pp.getarr("reference_velocity", m_ref_vel);

    // Which coordinate directions to force
    pp.queryarr("force_coord_directions", m_fcoord);

    // Based upon Allaerts & Meyers (JFM, 2017) and Durran & Klemp (AMS, 1983)
    pp.query("lateral_damping_start", m_meso_start);
    pp.query("west_damping_length", m_west_damping_len);
    pp.query("east_damping_length", m_east_damping_len);
    pp.query("north_damping_length", m_north_damping_len);
    pp.query("south_damping_length", m_south_damping_len);
    pp.query("vertical_cutoff", m_vertical_cutoff);
}

RayleighDamping::~RayleighDamping() = default;

void RayleighDamping::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState fstate,
    const amrex::Array4<amrex::Real>& src_term) const
{

    const auto& problo = m_mesh.Geom(lev).ProbLoArray();
    const auto& probhi = m_mesh.Geom(lev).ProbHiArray();
    const auto& dx = m_mesh.Geom(lev).CellSizeArray();

    const amrex::Real tau = m_tau;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> ref_vel{
        m_ref_vel[0], m_ref_vel[1], m_ref_vel[2]};
    const auto& vel =
        m_velocity.state(field_impl::dof_state(fstate))(lev).const_array(mfi);

    // Constants used to determine the fringe region coefficient
    const amrex::Real dRD = m_dRD;
    const amrex::Real dFull = m_dFull;

    // Which coordinate directions to force
    const amrex::Real fx = m_fcoord[0];
    const amrex::Real fy = m_fcoord[1];
    const amrex::Real fz = m_fcoord[2];
    //! Lateral
    const amrex::Real meso_sponge_start = m_meso_start;
    const amrex::Real west_damping_len = m_west_damping_len;
    const amrex::Real east_damping_len = m_east_damping_len;
    const amrex::Real north_damping_len = m_north_damping_len;
    const amrex::Real south_damping_len = m_south_damping_len;
    const amrex::Real vertical_cutoff = m_vertical_cutoff;
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        amrex::Real coeff = 0.0;
        const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
        const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
        const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

        if (probhi[2] - z > dRD + dFull) {
            coeff = 0.0;
        } else if (probhi[2] - z > dFull) {
            coeff = 0.5 * std::cos(M_PI * (probhi[2] - dFull - z) / dRD) + 0.5;
        } else {
            coeff = 1.0;
        }
        const amrex::Real west_zrd = problo[0] + west_damping_len;
        const amrex::Real west_Rayleigh_coeff =
            std::sin(0.5 * M_PI * (x - west_zrd) / (problo[0] - west_zrd));
        const amrex::Real damp_west =
            (x < west_zrd) ? west_Rayleigh_coeff * west_Rayleigh_coeff : 0.0;
        const amrex::Real east_zrd = probhi[0] - east_damping_len;
        const amrex::Real east_Rayleigh_coeff =
            std::sin(0.5 * M_PI * (x - east_zrd) / (probhi[0] - east_zrd));
        const amrex::Real damp_east =
            (x > east_zrd) ? east_Rayleigh_coeff * east_Rayleigh_coeff : 0.0;
        const amrex::Real south_zrd = problo[1] + south_damping_len;
        const amrex::Real south_Rayleigh_coeff =
            std::sin(0.5 * M_PI * (y - south_zrd) / (problo[1] - south_zrd));
        const amrex::Real damp_south =
            (y < south_zrd) ? south_Rayleigh_coeff * south_Rayleigh_coeff : 0.0;
        const amrex::Real north_zrd = probhi[1] - north_damping_len;
        const amrex::Real north_Rayleigh_coeff =
            std::sin(0.5 * M_PI * (y - north_zrd) / (probhi[1] - north_zrd));
        const amrex::Real damp_north =
            (y > north_zrd) ? north_Rayleigh_coeff * north_Rayleigh_coeff : 0.0;
        const amrex::Real zhigh = vertical_cutoff;
        const amrex::Real Rayleigh_coeff = std::sin(
            0.5 * M_PI * (z - meso_sponge_start) / (zhigh - meso_sponge_start));
        const amrex::Real horizontal_Rayleigh_z =
            (z > meso_sponge_start && z < zhigh)
                ? Rayleigh_coeff * Rayleigh_coeff *
                      (ref_vel[2] - vel(i, j, k, 2)) / tau
                : 0.0;
        src_term(i, j, k, 0) +=
            fx * coeff * (ref_vel[0] - vel(i, j, k, 0)) / tau;
        src_term(i, j, k, 1) +=
            fy * coeff * (ref_vel[1] - vel(i, j, k, 1)) / tau;
        src_term(i, j, k, 2) +=
            fz * coeff * (ref_vel[2] - vel(i, j, k, 2)) / tau +
            std::max(damp_west + damp_east + damp_north + damp_south, 1.0) *
                horizontal_Rayleigh_z;
    });
}

} // namespace amr_wind::pde::icns

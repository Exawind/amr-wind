#include "amr-wind/equation_systems/icns/source_terms/IdealizedHurricaneForcing.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/core/vs/vstraits.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"

namespace amr_wind::pde::icns {

IdealizedHurricaneForcing::IdealizedHurricaneForcing(const CFDSim& sim)
    : m_velocity(sim.repo().get_field("velocity"))
{

    amrex::Real coriolis_factor;
    {
        // Read the rotational time period (in seconds)
        amrex::ParmParse pp("CoriolisForcing");
        amrex::Real rot_time_period = 86400.0;
        pp.query("rotational_time_period", rot_time_period);
        coriolis_factor = 2.0 * utils::two_pi() / rot_time_period;
        amrex::Print() << "Geostrophic forcing: Coriolis factor = "
                       << coriolis_factor << std::endl;
        amrex::Real latitude = 90.0;
        pp.query("latitude", latitude);
        AMREX_ALWAYS_ASSERT(
            amrex::Math::abs(latitude - 90.0) <
            static_cast<amrex::Real>(vs::DTraits<float>::eps()));
    }

    {
        // Read the geostrophic wind speed vector (in m/s)
        amrex::ParmParse pp("IdealizedHurricaneForcing");
        pp.query("gradient_wind", m_Ug);
        pp.query("eyewall_radial_distance", m_R);
    }
}

IdealizedHurricaneForcing::~IdealizedHurricaneForcing() = default;

void IdealizedHurricaneForcing::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState fstate,
    const amrex::Array4<amrex::Real>& src_term) const
{

    const amrex::Real R = m_R;
    const amrex::Real Ug = m_Ug;

    const auto& vel =
        m_velocity.state(field_impl::dof_state(fstate))(lev).const_array(mfi);

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const amrex::Real umean = 0.0;
        const amrex::Real dumeandR = 1.0;
        const amrex::Real vmean = 0.0;

        src_term(i, j, k, 0) +=
            vel(i, j, k, 1) * umean / R + vel(i, j, k, 0) * dumeandR;
        src_term(i, j, k, 1) += -vel(i, j, k, 1) * vmean / R -
                                umean * vel(i, j, k, 0) / R + Ug * Ug / R;
        src_term(i, j, k, 2) += 0.0;
    });
}

} // namespace amr_wind::pde::icns

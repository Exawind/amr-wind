#include "amr-wind/equation_systems/icns/source_terms/BoussinesqBuoyancy.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/FieldUtils.H"

#include "AMReX_ParmParse.H"

namespace amr_wind::pde::icns {

/** Boussinesq buoyancy source term for ABL simulations
 *
 *  Reads in the following parameters from `BoussinesqBuoyancy` namespace:
 *
 *  - `reference_temperature` (Mandatory) temperature (`T0`) in Kelvin
 *  - `thermal_expansion_coeff` Optional, default = `1.0 / T0`
 *  - `gravity` acceleration due to gravity (m/s)
 */
BoussinesqBuoyancy::BoussinesqBuoyancy(const CFDSim& sim)
    : m_temperature(sim.repo().get_field("temperature"))
{
    amrex::ParmParse pp_boussinesq_buoyancy(
        amr_wind::pde::icns::BoussinesqBuoyancy::identifier());
    pp_boussinesq_buoyancy.get("reference_temperature", m_ref_theta);

    std::string transport_model_name = "ConstTransport";
    {
        amrex::ParmParse pp("transport");
        pp.query("model", transport_model_name);
    }
    m_transport = transport::TransportModel::create(transport_model_name, sim);

    m_is_vof = sim.repo().field_exists("vof");
    if (m_is_vof) {
        m_vof = &sim.repo().get_field("vof");
    } else {
        // Point to something, will not be used
        m_vof = &m_temperature;
    }

    // gravity in `incflo` namespace
    amrex::ParmParse pp_incflo("incflo");
    pp_incflo.queryarr("gravity", m_gravity);
}

BoussinesqBuoyancy::~BoussinesqBuoyancy() = default;

void BoussinesqBuoyancy::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState fstate,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const amrex::Real T0 = m_ref_theta;
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{
        m_gravity[0], m_gravity[1], m_gravity[2]};

    const bool ivf = m_is_vof;
    const auto& vof_arr = (*m_vof)(lev).const_array(mfi);
    constexpr amrex::Real tol = 1e-12;

    const auto& temp =
        m_temperature.state(field_impl::phi_state(fstate))(lev).const_array(
            mfi);
    const auto& beta = (*m_transport->beta())(lev).const_array(mfi);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const amrex::Real T = temp(i, j, k, 0);
        const amrex::Real fac_air = beta(i,j,k) * (T0 - T);
        // If vof exists, ignore Boussinesq term in cells with liquid
        // If no vof, assume single phase and use the term for air everywhere
        const amrex::Real fac =
            ivf ? (vof_arr(i, j, k) > tol ? 0.0 : fac_air) : fac_air;

        src_term(i, j, k, 0) += gravity[0] * fac;
        src_term(i, j, k, 1) += gravity[1] * fac;
        src_term(i, j, k, 2) += gravity[2] * fac;
    });
}

} // namespace amr_wind::pde::icns

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
    std::string transport_model_name = "ConstTransport";
    {
        amrex::ParmParse pp("transport");
        pp.query("model", transport_model_name);
    }
    m_transport = transport::TransportModel::create(transport_model_name, sim);
    m_beta = m_transport->beta();

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
    const amrex::Real T0 = m_transport->reference_temperature();
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{
        m_gravity[0], m_gravity[1], m_gravity[2]};

    const auto& temp =
        m_temperature.state(field_impl::phi_state(fstate))(lev).const_array(
            mfi);
    const auto& beta = (*m_beta)(lev).const_array(mfi);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const amrex::Real T = temp(i, j, k, 0);
        const amrex::Real fac = beta(i, j, k) * (T0 - T);

        src_term(i, j, k, 0) += gravity[0] * fac;
        src_term(i, j, k, 1) += gravity[1] * fac;
        src_term(i, j, k, 2) += gravity[2] * fac;
    });
}

} // namespace amr_wind::pde::icns

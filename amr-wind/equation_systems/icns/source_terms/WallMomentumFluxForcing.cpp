#include "amr-wind/equation_systems/icns/source_terms/WallMomentumFluxForcing.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/FieldUtils.H"
#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/wind_energy/ShearStress.H"

#include "AMReX_ParmParse.H"

namespace amr_wind::pde::icns {

// FIXME: comments out of date
/** Boussinesq buoyancy source term for ABL simulations
 *
 *  Reads in the following parameters from `ABLMeanBoussinesq` namespace:
 *
 *  - `reference_temperature` (Mandatory) temperature (`T0`) in Kelvin
 *  - `thermal_expansion_coeff` Optional, default = `1.0 / T0`
 *  - `gravity` acceleration due to gravity (m/s)
 *  - `read_temperature_profile`
 *  - `tprofile_filename`
 */
WallMomentumFluxForcing::WallMomentumFluxForcing(const CFDSim& sim)
    : m_mesh(sim.mesh())
    , m_velocity(sim.repo().get_field("velocity"))
    , m_density(sim.repo().get_field("density"))
    , m_mo(sim.physics_manager().get<amr_wind::ABL>().abl_wall_function().mo())
    , m_wall_momentum_flux_source(
          sim.repo().get_field("wall_shear_stress_src_term"))
{

    // some parm parse stuff?
    amrex::ParmParse pp("ABL");
    pp.query("wall_shear_stress_type", m_wall_shear_stress_type);
    pp.query("normal_direction", m_direction);
    AMREX_ASSERT((0 <= m_direction) && (m_direction < AMREX_SPACEDIM));
}

WallMomentumFluxForcing::~WallMomentumFluxForcing() = default;

template <typename ShearStress>
void WallMomentumFluxForcing::compute_wall_src(
    const ShearStress& tau,
    const amrex::Real& weightLow,
    const amrex::Real& weightHigh,
    const amrex::Array4<amrex::Real>& src_term,
    const amrex::Array4<amrex::Real>& plotField,
    const amrex::Array4<amrex::Real>& velocityField)
{
    amrex::ParallelFor(amrex::bdryLo(bx, idir),
    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        // Get the local velocity at the cell center adjacent
        // to this wall face.
        const amrex::Real uLow = velocityField(i, j, kLow, 0);
        const amrex::Real vLow = velocityField(i, j, kLow, 1);
        
        const amrex::Real uHigh = velocityField(i, j, kHigh, 0);
        const amrex::Real vHigh = velocityField(i, j, kHigh, 1);
        
        const amrex::Real u = weightLow * uLow + weightHigh * uHigh;
        const amrex::Real v = weightLow * vLow + weightHigh * vHigh;
        const amrex::Real S = std::sqrt((u * u) + (v * v));
        
        // Get local tau_wall based on the local conditions and
        // mean state based on Monin-Obukhov similarity.
        amrex::Real tau_xz = tau.calc_vel_x(u, S);
        amrex::Real tau_yz = tau.calc_vel_y(v, S);
        
        // Adding the source term as surface stress vector times surface
        // area divided by cell volume (division by cell volume is to make
        // this a source per unit volume).
        plotField(i, j, k, 0) = -tau_xz;
        plotField(i, j, k, 1) = -tau_yz;
        plotField(i, j, k, 2) = 0.0;
        
        src_term(i, j, k, 0) -= (tau_xz * dx[0] * dx[1]) / dV;
        src_term(i, j, k, 1) -= (tau_yz * dx[0] * dx[1]) / dV;
        src_term(i, j, k, 2) += 0.0;
    });
}

void WallMomentumFluxForcing::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState fstate,
    const amrex::Array4<amrex::Real>& src_term) const
{
    // Overall geometry information
    const auto& geom = m_mesh.Geom(lev);

    // Mesh cell size information
    const auto& dx = m_mesh.Geom(lev).CellSizeArray();
    amrex::Real dV = 1.0;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        dV *= dx[dir];
    }

    // Domain size information.
    const auto& domain = geom.Domain();

    //
    const int idir = m_direction;

    const auto& velocityField =
        m_velocity.state(field_impl::dof_state(fstate))(lev).const_array(mfi);

    auto plotField = m_wall_momentum_flux_source(lev).array(mfi);

    // FieldState densityState = field_impl::phi_state(fstate);
    // const auto& density =
    // m_density.state(densityState)(lev).const_array(mfi);

    // Get the desired sampling height.
    const amrex::Real zref = m_mo.zref;

    // Figure out on what grid level the sampling should occur and the
    // interpolation weights.
    const amrex::Real index = (zref / dx[idir]) - 0.5;

    int kLow = int(std::floor(index));
    int kHigh = int(std::ceil(index));
    // if kLow and kHigh point to the same grid cell, separate them by one cell.
    if (kLow == kHigh) {
        kHigh += 1;
    }
    // if kLow lies below the grid box, bump it up to the first grid level.
    if (kLow < 0) {
        kLow = 0;
        kHigh = kLow + 1;
    }
    // if kHigh lies outside the grid box, bump it down to the last grid level.
    else if (kHigh > bx.bigEnd(idir)) {
        kHigh = bx.bigEnd(idir);
        kLow = kHigh - 1;
    }
    // kLow = (kLow < 0) ? 0 : kLow;
    // kHigh = (kLow == kHigh) ? kHigh + 1 : kHigh;
    // kHigh = (kHigh > bx.bigEnd(idir)) ? bx.bigEnd(idir) : kHigh;
    // kLow = (kLow == kHigh) ? kLow - 1 : kLow;

    const amrex::Real weightLow = amrex::Real(kHigh) - index;
    const amrex::Real weightHigh = index - amrex::Real(kLow);

    /*
        std::cout << "index = " << index << std::endl;
        std::cout << "kLow = " << kLow << std::endl;
        std::cout << "kHigh = " << kHigh << std::endl;
        std::cout << "weightLow = " << weightLow << std::endl;
        std::cout << "weightHigh = " << weightHigh << std::endl;
        std::cout << "bx.smallEnd(0) : bx.bigEnd(0) = " << bx.smallEnd(0) << " :
       " << bx.bigEnd(0) << std::endl; std::cout << "bx.smallEnd(1) :
       bx.bigEnd(1) = " << bx.smallEnd(1) << " : " << bx.bigEnd(1) << std::endl;
        std::cout << "bx.smallEnd(2) : bx.bigEnd(2) = " << bx.smallEnd(2) << " :
       " << bx.bigEnd(2) << std::endl;
    */

    if (!(bx.smallEnd(idir) == domain.smallEnd(idir))) return;
    if (idir != 2) return;

    if (m_wall_shear_stress_type == "constant") {
        auto tau = ShearStressConstant(m_mo);
        compute_wall_src(tau, weightLow, weightHigh, src_term, plotField, velocityField);
    } else if (m_wall_shear_stress_type == "default") {
        auto tau = ShearStressDefault(m_mo);
        compute_wall_src(tau, weightLow, weightHigh, src_term, plotField, velocityField);
    } else if (m_wall_shear_stress_type == "local") {
        auto tau = ShearStressLocal(m_mo);
        compute_wall_src(tau, weightLow, weightHigh, src_term, plotField, velocityField);
    } else if (m_wall_shear_stress_type == "schumann") {
        auto tau = ShearStressSchumann(m_mo);
        compute_wall_src(tau, weightLow, weightHigh, src_term, plotField, velocityField);
    } else {
        auto tau = ShearStressMoeng(m_mo);
        compute_wall_src(tau, weightLow, weightHigh, src_term, plotField, velocityField);
    }
}
} // namespace amr_wind::pde::icns

#ifndef WALLMOMENTUMFLUXFORCING_H
#define WALLMOMENTUMFLUXFORCING_H

#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/wind_energy/MOData.H"

namespace amr_wind::pde::icns {

// class Field;

// FIXME: outdated comments
/** A forcing term that acts on the layer of cells adjacent to the wall
 *  to apply a desired wall shear stress.
 *  \ingroup icns_src we_abl
 *
 *  \f[
 *    \vec{S} = \vec{tau_{wall}} * dx * dy
 *  \f]

 *  where tau_wall comes from a model based on Monin-Obukhov similarity
 *  theory and the definition of drag as \vec{D} = C_d * |\vec{V}| * \vec{V}
 */
class WallMomentumFluxForcing
    : public MomentumSource::Register<WallMomentumFluxForcing>
{
public:
    static const std::string identifier() { return "WallMomentumFluxForcing"; }

    explicit WallMomentumFluxForcing(const CFDSim& sim);

    ~WallMomentumFluxForcing() override;

    template <typename ShearStress>
    void compute_wall_src(
       const int& idir,
       const int& kLow,
       const int& kHigh,
       const amrex::Real& dV,
       const amrex::Real& weightLow,
       const amrex::Real& weightHigh,
       const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>& dx,
       const amrex::Box& bx,
       const ShearStress& tau,
       const amrex::Array4<amrex::Real>& src_term,
       const amrex::Array4<amrex::Real>& plotField,
       const amrex::Array4<const amrex::Real>& velocityField) const;
        
    void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const;

private:
    const amrex::AmrCore& m_mesh;
    const Field& m_velocity;
    const Field& m_density;

    //! Reference to Monin-Obukhov
    const MOData& m_mo;

    //! Wall-normal direction
    int m_direction{2};

    //! Particular shear stress model
    std::string m_wall_shear_stress_type{"moeng"};

    Field& m_wall_momentum_flux_source;
};
} // namespace amr_wind::pde::icns

#endif /* WALLMOMENTUMFLUXFORCING_H */

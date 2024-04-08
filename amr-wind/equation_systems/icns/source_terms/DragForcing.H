#ifndef DRAGFORCING_H
#define DRAGFORCING_H

#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind::pde::icns {

/** Adding Terrain-induced drag forcing to momentum equation ..
 *
 *  \ingroup icns_src DragForcing
 *
 *  \sa DragForcing
 */
class DragForcing : public MomentumSource::Register<DragForcing>
{
public:
    static std::string identifier() { return "DragForcing"; }

    explicit DragForcing(const CFDSim& sim);

    ~DragForcing() override;

    void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    const CFDSim& m_sim;
    const amrex::AmrCore& m_mesh;
    const Field& m_velocity;
    Field& m_terrainBlank;
    Field& m_terrainDrag;
    // amrex::Gpu::DeviceVector<amrex::Real> device_vel_ht;
    // amrex::Gpu::DeviceVector<amrex::Real> device_vel_vals;
    amrex::Real m_dragCoeff{8.0};
    amrex::Real m_spongeStrength{1.0};
    amrex::Real m_spongeDensity{1.0};
    amrex::Real m_spongePercentX{20};
    amrex::Real m_spongePercentY{20};
};

} // namespace amr_wind::pde::icns

#endif /* DRAGFORCING_H */

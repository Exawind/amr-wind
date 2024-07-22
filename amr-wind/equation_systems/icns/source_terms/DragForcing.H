#ifndef DRAGFORCING_H
#define DRAGFORCING_H

#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind::pde::icns {

/** Adds the forcing term to include the presence of immersed boundary
 *
 *  \ingroup icns_src
 *
 *
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
    const SimTime& m_time;
    const CFDSim& m_sim;
    const amrex::AmrCore& m_mesh;
    const Field& m_velocity;
    amrex::Gpu::DeviceVector<amrex::Real> m_device_vel_ht;
    amrex::Gpu::DeviceVector<amrex::Real> m_device_vel_vals;
    // Smaller value gives a nice smooth start 
    amrex::Real m_drag{10.0};
    amrex::Real m_sponge_strength{1.0};
    amrex::Real m_sponge_density{1.0};
    amrex::Real m_sponge_distanceX{1000};
    amrex::Real m_sponge_distanceY{1000};
};

} // namespace amr_wind::pde::icns

#endif

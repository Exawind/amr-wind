#ifndef MULTIPHASEDRAGFORCING_H
#define MULTIPHASEDRAGFORCING_H

#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind::pde::icns {

/** Adds the forcing term to include the presence of immersed boundary for
 * multiphase flows
 *
 *  \ingroup icns_src
 *
 *
 */
class MultiphaseDragForcing
    : public MomentumSource::Register<MultiphaseDragForcing>
{
public:
    static std::string identifier() { return "MultiphaseDragForcing"; }

    explicit MultiphaseDragForcing(const CFDSim& sim);

    ~MultiphaseDragForcing() override;

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
    const Field* m_target_vel{nullptr};
    const Field* m_target_levelset{nullptr};
};

} // namespace amr_wind::pde::icns

#endif

#ifndef GEOSTROPHICFORCING_H
#define GEOSTROPHICFORCING_H

#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/core/SimTime.H"

namespace amr_wind::pde::icns {

/** Apply source term equivalent to drive the flow to match given geostrophic
 * wind \ingroup icns_src we_abl
 */
class GeostrophicForcing : public MomentumSource::Register<GeostrophicForcing>
{
public:
    static std::string identifier() { return "GeostrophicForcing"; }

    explicit GeostrophicForcing(const CFDSim& /*unused*/);

    ~GeostrophicForcing() override;

    void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    //! Time
    const SimTime& m_time;

    const amrex::AmrCore& m_mesh;

    //! Activated when water is present in domain
    bool m_use_phase_ramp{false};

    //! Number of cells in band to prevent forcing near liquid
    int m_n_band{2};

    //! Target velocity
    amrex::Vector<amrex::Real> m_target_vel{{0.0, 0.0, 0.0}};

    //! Forcing source term (pressure gradient)
    amrex::Vector<amrex::Real> m_g_forcing{{0.0, 0.0, 0.0}};

    bool m_is_horizontal{false};

    //! Bool for writing forcing time table
    bool m_write_force_timetable{false};
    //! File name for forcing time table output
    std::string m_force_timetable;
    //! Output frequency for forces
    int m_force_outfreq{1};
    //! Output start time for force
    amrex::Real m_force_outstart{0.0};

    //! Height from water interface where force is off
    amrex::Real m_forcing_mphase0;
    //! Height from water interface over which force is ramped
    amrex::Real m_forcing_mphase1;

    //! Local storage of interface location
    amrex::Real m_water_level;

    //! VOF field, to avoid forcing on liquid above force-off height
    const Field* m_vof;
};

} // namespace amr_wind::pde::icns

#endif /* GEOSTROPHICFORCING_H */

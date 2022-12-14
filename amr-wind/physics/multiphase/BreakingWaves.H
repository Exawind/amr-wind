#ifndef NWB_H
#define NWB_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/SimTime.H"

/**
 *  \defgroup Breaking Waves modeling
 *  Breaking Waves modeling
 *
 *  This group documents all the wave mechanics related physics models available
 * in AMR-Wind
 *
 *  \ingroup physics
 */

/**
 *  \defgroup we_nwb NWB
 *  Numerical Wave Basin (NWB) modeling
 *
 *  \ingroup waves
 */

namespace amr_wind {

/** Numerical Wave Basin physics
 *  \ingroup we_abl
 *
 */
class BreakingWaves : public Physics::Register<BreakingWaves>
{
public:
    static std::string identifier() { return "BreakingWaves"; }

    explicit BreakingWaves(CFDSim& sim);

    ~BreakingWaves() override;

    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void post_init_actions() override {}

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override {}

private:
    const CFDSim& m_sim;

    Field& m_velocity;
    Field& m_levelset;
    Field& m_density;

    //! initial tke value
    amrex::Real m_tke0{0.000001};

    //! initial sdr value
    amrex::Real m_sdr0{1.0};

    //! Initial free surface amplitude magnitude
    amrex::Real m_amplitude{0.1};

    amrex::Real m_air_vel_shear_rate{0.1};

    //! Initial free surface wavelength
    amrex::Real m_wavelength{2.0};

    //! Initial zero-level free-surface water depth
    amrex::Real m_waterlevel{0.0};

    bool m_laminar{false};
};

} // namespace amr_wind

#endif /* NWB_H */

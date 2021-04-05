#ifndef NWB_H
#define NWB_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/SimTime.H"

/**
 *  \defgroup Water Waves modeling
 *  Water Waves modeling
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
class NWB : public Physics::Register<NWB>
{
public:
    static const std::string identifier() { return "NWB"; }

    explicit NWB(CFDSim& sim);

    virtual ~NWB();

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

    //! Initial free surface amplitude magnitude
    amrex::Real m_amplitude{0.1};

    //! Initial free surface wavelength
    amrex::Real m_wavelength{2.0};

    //! Initial zero-level free-surface water depth
    amrex::Real m_waterlevel{0.0};

    //! Airflow velocity magnitude
    amrex::Real m_airflow_velocity{1.0};
};

} // namespace amr_wind

#endif /* NWB_H */

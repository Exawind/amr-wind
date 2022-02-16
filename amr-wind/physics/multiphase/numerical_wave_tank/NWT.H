#ifndef NWT_H
#define NWT_H

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
 *  \defgroup we_nwt NWT
 *  Numerical Wave Tank (NWT) modeling
 *
 *  \ingroup waves
 */

namespace amr_wind {

/** Numerical Wave Tank physics
 *  \ingroup we_nwt
 *
 */
class NWT : public Physics::Register<NWT>
{
public:
    static const std::string identifier() { return "NWT"; }

    explicit NWT(CFDSim& sim);

    virtual ~NWT();

    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void post_init_actions() override {}

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override;

    void apply_relaxation_method();

private:
    const CFDSim& m_sim;

    Field& m_velocity;
    Field& m_levelset;
    Field& m_density;
    Field& m_vof;

    // bool m_has_ramp;

    // Enumerator for wave generation theory
    enum class wave_generator {
        LinearWaves,
        StokesWaves,
        BiChromatic,
        Spectrum
    };

    //! Initial free surface amplitude magnitude
    amrex::Real m_amplitude{0.1};

    //! Initial free surface wavelength
    amrex::Real m_wavelength{1.0};

    //! Initial zero-sea-level
    amrex::Real m_zsl{0.0};

    //! Water depth
    amrex::Real m_waterdepth{1.0};

    //! Wavemaker ramp function time-period
    amrex::Real m_ramp_period{10.0};

    //! Generation length
    amrex::Real m_gen_length{1.0};

    //! Numerical beach absorption length
    amrex::Real m_beach_length{2.0};

    //! Numerical beach absorption length factor
    amrex::Real m_absorb_length_factor{1.0};

    //! Numerical beach start location
    amrex::Real m_x_start_beach{10};
};

} // namespace amr_wind

#endif /* NWT_H */

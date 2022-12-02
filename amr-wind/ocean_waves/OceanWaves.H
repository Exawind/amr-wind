#ifndef OCEAN_WAVES_H
#define OCEAN_WAVES_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/FieldRepo.H"

/** \defgroup Ocean Waves module
 *  Ocean Waves class
 *  \ingroup ocean_waves
 */

namespace amr_wind {

class Field;

namespace ocean_waves {

class OceanWavesModel;

/** ocean waves modeling
 *
 *  \ingroup ocean waves
 *
 *  This class provides an interface to model
 *
 *  \sa WaveModel
 */
class OceanWaves : public Physics::Register<OceanWaves>
{
public:
    static std::string identifier() { return "OceanWaves"; }

    explicit OceanWaves(CFDSim& /*sim*/);

    ~OceanWaves() override;

    void
    initialize_fields(int /*level*/, const amrex::Geometry& /*geom*/) override;

    void pre_init_actions() override;

    void post_init_actions() override;

    void post_regrid_actions() override;

    void pre_advance_work() override;

    void post_advance_work() override;

protected:
    virtual void prepare_outputs();

private:
    void relaxation_zones();

    CFDSim& m_sim;

    //! Unique pointer to the ocean waves model
    std::unique_ptr<OceanWavesModel> m_owm;

    //! Ocean waves target free surface levelset function
    Field& m_ow_levelset;

    //! Ocean waves target volume-of-fluid
    Field& m_ow_vof;

    //! Ocean waves target velocity
    Field& m_ow_velocity;

    //! Ocean waves pressure gradient
    // Field& m_ow_pressure;
};

} // namespace ocean_waves
} // namespace amr_wind

#endif /* OCEAN_WAVES_H */

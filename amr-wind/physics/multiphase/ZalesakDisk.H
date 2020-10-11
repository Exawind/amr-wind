#ifndef ZalesakDisk_H
#define ZalesakDisk_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/physics/multiphase/ZalesakDiskFieldInit.H"
#include "amr-wind/core/Field.H"

/** Multiphase Zalesak disk physics
 *  \ingroup mulitphase_physics
 *
 */

namespace amr_wind {

class ZalesakDisk : public Physics::Register<ZalesakDisk>
{
public:
    static const std::string identifier() { return "ZalesakDisk"; }

    ZalesakDisk(CFDSim& sim);

    virtual ~ZalesakDisk() = default;

    //! Initialize the levelset and velocity fields for Vortex Patch simulations

    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void post_init_actions() override {}

    void post_regrid_actions() override {}

    void pre_advance_work() override;

    void post_advance_work() override;

private:
    const CFDSim& m_sim;

    Field& m_velocity;
    Field& m_levelset;

    //! ZalesakDisk field initializer instance
    std::unique_ptr<ZalesakDiskFieldInit> m_field_init;

    //! prescribed velocity period
    amrex::Real m_TT{6.0};
};

} // namespace amr_wind

#endif /* ZalesakDisk_H */

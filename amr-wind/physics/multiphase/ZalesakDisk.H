#ifndef ZalesakDisk_H
#define ZalesakDisk_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"

/** Multiphase Zalesak disk physics
 *  \ingroup mulitphase_physics
 *
 */

namespace amr_wind {

class ZalesakDisk : public Physics::Register<ZalesakDisk>
{

    static_assert(
        AMREX_SPACEDIM == 3, "ZalesakDisk requires 3 dimensional mesh");

public:
    static const std::string identifier() { return "ZalesakDisk"; }

    explicit ZalesakDisk(CFDSim& sim);

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
    Field& m_density;

    //! Initial ZalesakDisk location
    amrex::Vector<amrex::Real> m_loc{{0.5, 0.72, 0.24}};

    //! sphere radius value
    amrex::Real m_radius{0.16};

    //! slot width value
    amrex::Real m_width{0.04};

    //! slot depth
    amrex::Real m_depth{0.2};

    //! prescribed velocity period
    amrex::Real m_TT{6.0};
};

} // namespace amr_wind

#endif /* ZalesakDisk_H */

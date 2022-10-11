#ifndef ZalesakDiskScalarVel_H
#define ZalesakDiskScalarVel_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"

/** Multiphase Zalesak disk physics
 *  \ingroup mulitphase_physics
 *
 */
namespace amr_wind::zds {

namespace {
struct SCexact
{
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real operator()(
        amrex::Real /*xc*/,
        amrex::Real /*yc*/,
        amrex::Real /*zc*/,
        amrex::Real /*x*/,
        amrex::Real /*y*/,
        amrex::Real /*z*/,
        amrex::Real /*radius*/,
        amrex::Real /*alpha*/) const;
    const int m_comp{0};
};
} // namespace

class ZalesakDiskScalarVel : public Physics::Register<ZalesakDiskScalarVel>
{

    static_assert(
        AMREX_SPACEDIM == 3,
        "ZalesakDiskScalarVel requires 3 dimensional mesh");

public:
    static std::string identifier() { return "ZalesakDiskScalarVel"; }

    explicit ZalesakDiskScalarVel(CFDSim& sim);

    ~ZalesakDiskScalarVel() override = default;

    //! Initialize the levelset and velocity fields for Vortex Patch simulations

    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void post_init_actions() override {}

    void post_regrid_actions() override {}

    void pre_advance_work() override;

    void post_advance_work() override;

    template <typename T>
    amrex::Real compute_error(const Field& /*field*/);

private:
    const CFDSim& m_sim;

    Field& m_velocity;
    Field& m_levelset;
    Field& m_density;

    //! Initial ZalesakDiskScalarVel location
    amrex::Vector<amrex::Real> m_loc{{0.5, 0.72, 0.24}};

    //! sphere radius value
    amrex::Real m_radius{0.16};

    //! slot width value
    amrex::Real m_width{0.04};

    //! slot depth
    amrex::Real m_depth{0.2};

    //! prescribed velocity period
    amrex::Real m_TT{6.0};

    void output_error();

    //! output precision
    const int m_w{18};

    //! error log file
    std::string m_output_fname{"zds.log"};
};

} // namespace amr_wind::zds

#endif /* ZalesakDiskScalarVel_H */

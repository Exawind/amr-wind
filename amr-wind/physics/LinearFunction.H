#ifndef LinearFunction_H
#define LinearFunction_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"

namespace amr_wind {

/** Taylor Green Vortex physics
 *  \ingroup physics
 */
class LinearFunction : public Physics::Register<LinearFunction>
{
public:
    static const std::string identifier() { return "LinearFunction"; }

    explicit LinearFunction(const CFDSim& sim);

    virtual ~LinearFunction() = default;
    //! Initialize the temperature and velocity fields for BoussinesqBubble
    //! simulations
    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void post_init_actions() override {}

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override {}

private:
    const FieldRepo& m_repo;
    Field& m_velocity;

    //! initial density value
    amrex::Real m_rho{1.0};
    amrex::Real m_xcoef{0.0};
    amrex::Real m_ycoef{0.0};
    amrex::Real m_zcoef{0.0};
};

} // namespace amr_wind

#endif /* LinearFunction_H */

#ifndef ScalarAdvection_H
#define ScalarAdvection_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"

namespace amr_wind {

/** Boussinesq Bubble physics
 *  \ingroup physics
 *
 *  \sa ScalarAdvectionFieldInit
 */
class ScalarAdvection : public Physics::Register<ScalarAdvection>
{
public:
    static std::string identifier() { return "ScalarAdvection"; }

    explicit ScalarAdvection(CFDSim& sim);

    ~ScalarAdvection() override = default;

    //! Initialize the temperature and velocity fields for ScalarAdvection
    //! simulations
    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void post_init_actions() override {}

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override {}

private:
    Field& m_velocity;
    Field& m_density;

    //! Temperature field
    Field* m_temperature{nullptr};

    //! Initial bubble location
    amrex::Vector<amrex::Real> m_loc{{0.25, 0.25, 0.5}};

    //! Initial density field
    amrex::Real m_rho{1.0};

    //! tracer inner value
    amrex::Real m_tracer_inner{0.0};

    //! tracer outer value
    amrex::Real m_tracer_outer{0.01};

    //! tracer bubble radius value
    amrex::Real m_tracer_radius{0.1};
};

} // namespace amr_wind

#endif /* ScalarAdvection_H */

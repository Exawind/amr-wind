#ifndef FreeStream_H
#define FreeStream_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"

namespace amr_wind {

/** Freestream physics
 */
class FreeStream : public Physics::Register<FreeStream>
{
public:
    static const std::string identifier() { return "FreeStream"; }

    FreeStream(const CFDSim& sim);

    virtual ~FreeStream() = default;

    //! Initialize the temperature and velocity fields for BoussinesqBubble simulations
    void initialize_fields(
        int level,
        const amrex::Geometry& geom) override;

    void post_init_actions() override {}

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override {}

private:

    Field& m_velocity;
    Field& m_density;

    //! initial density value
    amrex::Real m_rho{1.0};
    amrex::Vector<amrex::Real> m_vel{{1.0, 0.0, 0.0}};

};

}

#endif /* FreeStream_H */

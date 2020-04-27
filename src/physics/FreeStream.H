#ifndef FreeStream_H
#define FreeStream_H

#include "Physics.H"
#include "Field.H"

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
    amrex::Real m_u{1.0};
    amrex::Real m_v{0.0};
    amrex::Real m_w{0.0};

};

}

#endif /* FreeStream_H */

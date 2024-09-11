#ifndef Intermittency_H
#define Intermittency_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"

namespace amr_wind {

/** Dummy transition model for overset exchange of the gamma variable 
 *  \ingroup physics
 */
class Intermittency : public Physics::Register<Intermittency>
{
public:
    static std::string identifier() { return "Intermittency"; }

    explicit Intermittency(const CFDSim& sim);

    ~Intermittency() override = default;
    //! Initialize the intermittency (gamma) field
    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void post_init_actions() override {}

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override {}

private:
    Field& m_intermittency;

    //! initial density value
    amrex::Real m_gamma{1.0};
};

} // namespace amr_wind

#endif /* TaylorGreenVortex_H */

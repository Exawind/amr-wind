#ifndef RayleighTaylor_H
#define RayleighTaylor_H

#include <memory>

#include "Physics.H"
#include "RayleighTaylorFieldInit.H"
#include "Field.H"

namespace amr_wind {

/** Rayleigh Taylor physics
 */
class RayleighTaylor : public Physics::Register<RayleighTaylor>
{
public:
    static const std::string identifier() { return "RayleighTaylor"; }

    RayleighTaylor(const CFDSim& sim);

    virtual ~RayleighTaylor() = default;

    //! Initialize the density and velocity fields for RayleighTaylor
    //! simulations
    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void post_init_actions() override {}

    void pre_advance_work() override {}

private:
    Field& m_velocity;
    Field& m_density;

    //! RayleighTaylor field initializer instance
    std::unique_ptr<RayleighTaylorFieldInit> m_field_init;
};

} // namespace amr_wind

#endif /* RayleighTaylor_H */

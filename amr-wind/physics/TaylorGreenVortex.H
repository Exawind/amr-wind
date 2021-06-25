#ifndef TaylorGreenVortex_H
#define TaylorGreenVortex_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"

namespace amr_wind {

/** Taylor Green Vortex physics
 *  \ingroup physics
 */
class TaylorGreenVortex : public Physics::Register<TaylorGreenVortex>
{
public:
    static const std::string identifier() { return "TaylorGreenVortex"; }

    explicit TaylorGreenVortex(const CFDSim& sim);

    virtual ~TaylorGreenVortex() = default;
    //! Initialize the temperature and velocity fields for BoussinesqBubble
    //! simulations
    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void post_init_actions() override {}

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override {}

private:
    Field& m_velocity;
    Field& m_density;
    Field& m_mesh_fac;

    //! probhi in unmapped mesh
    amrex::Vector<amrex::Real> m_probhi_unmapped{{0.0, 0.0, 0.0}};

    //! initial density value
    amrex::Real m_rho{1.0};
};

} // namespace amr_wind

#endif /* TaylorGreenVortex_H */

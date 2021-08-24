#ifndef FatCoredVortexRing_H
#define FatCoredVortexRing_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"

namespace amr_wind {

/** Vortex ring physics
 *  \ingroup physics
 */
class FatCoredVortexRing : public Physics::Register<FatCoredVortexRing>
{
public:
    static const std::string identifier() { return "FatCoredVortexRing"; }

    explicit FatCoredVortexRing(const CFDSim& sim);

    virtual ~FatCoredVortexRing() = default;
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
    Field& m_density;

    //! initial density value
    amrex::Real m_rho{1.0};

    amrex::Real m_R{1.0};
    amrex::Real m_Gamma{1.0};
    amrex::Real m_alpha{0.54857674};
};

struct VorticityTheta
{
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real operator()(
        const amrex::Real,
        const amrex::Real,
        const amrex::Real,
        const amrex::Real,
        const amrex::Real) const;
    const int m_comp{1};
};

} // namespace amr_wind

#endif /* FatCoredVortexRing_H */

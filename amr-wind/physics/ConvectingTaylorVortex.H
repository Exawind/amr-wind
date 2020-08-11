#ifndef ConvectingTaylorVortex_H
#define ConvectingTaylorVortex_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind {
namespace ctv {

namespace {
struct UExact
{
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real operator()(
        const amrex::Real,
        const amrex::Real,
        const amrex::Real,
        const amrex::Real,
        const amrex::Real,
        const amrex::Real) const;
    const int m_comp{0};
};

struct VExact
{
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real operator()(
        const amrex::Real,
        const amrex::Real,
        const amrex::Real,
        const amrex::Real,
        const amrex::Real,
        const amrex::Real) const;
    const int m_comp{1};
};

} // namespace

/** Convecting Taylor Vortex physics
 *  \ingroup physics
 */
class ConvectingTaylorVortex : public Physics::Register<ConvectingTaylorVortex>
{
public:
    static const std::string identifier() { return "ConvectingTaylorVortex"; }

    ConvectingTaylorVortex(const CFDSim& sim);

    virtual ~ConvectingTaylorVortex() = default;

    void initialize_fields(int level, const amrex::Geometry& geom) override;

    template <typename T>
    amrex::Real compute_error(const Field&);

    void post_init_actions() override;

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override;

private:
    const amr_wind::SimTime& m_time;
    const FieldRepo& m_repo;
    const amrex::AmrCore& m_mesh;
    Field& m_velocity;
    Field& m_density;
    void output_error();

    //! initial density value
    amrex::Real m_rho{1.0};

    //! mean x-velocity value
    amrex::Real m_u0{1.0};

    //! mean y-velocity value
    amrex::Real m_v0{1.0};

    //! damping
    amrex::Real m_omega;

    //! output precision
    const int m_w{18};

    //! error log file
    const std::string m_output_fname{"ctv.log"};
};
} // namespace ctv
} // namespace amr_wind

#endif /* ConvectingTaylorVortex_H */

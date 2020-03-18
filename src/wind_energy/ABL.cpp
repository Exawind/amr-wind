#include "ABL.H"
#include "ABLFieldInit.H"
#include "incflo.H"

#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"
#include "PlaneAveraging.H"

namespace amr_wind {

ABL::ABL(const SimTime& time, incflo* incflo_in)
    : m_time(time), m_incflo(incflo_in)
{
    amrex::ParmParse pp("abl");

    pp.query("use_boussinesq", m_has_boussinesq);
    pp.query("coriolis_effect", m_has_coriolis);
    pp.query("abl_forcing", m_has_driving_dpdx);

    // Instantiate the ABL field initializer
    m_field_init = std::make_unique<ABLFieldInit>();
    m_abl_wall_func = std::make_unique<ABLWallFunction>();

    if (m_has_driving_dpdx)
        m_abl_forcing = std::make_unique<ABLForcing>(m_time);
}

void ABL::initialize_fields(
    const amrex::Geometry& geom,
    LevelData& leveldata) const
{
    auto& velocity = leveldata.velocity;
    auto& density = leveldata.density;
    auto& scalars = leveldata.tracer;

    for (amrex::MFIter mfi(density); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();

        (*m_field_init)(
            vbx, geom, velocity.array(mfi), density.array(mfi),
            scalars.array(mfi));
    }
}

void ABL::add_momentum_sources(
    const amrex::Geometry& geom,
    const LevelData& leveldata,
    amrex::MultiFab& vel_forces) const
{
}

void ABL::pre_advance_work()
{
    // Spatial averaging on z-planes
    constexpr int direction = 2;
    const auto& geom = m_incflo->Geom(0);

    // TODO: Promote this to a class member
    PlaneAveraging pa(geom, *m_incflo->leveldata_vec()[0], direction);
    {
        // First cell height
        const amrex::Real fch = geom.ProbLo(direction) + 0.5 * geom.CellSize(direction);
        const amrex::Real vx = pa.line_velocity_xdir(fch);
        const amrex::Real vy = pa.line_velocity_ydir(fch);

        const amrex::Real uground = std::sqrt(vx * vx + vy * vy);
        const amrex::Real utau = m_abl_wall_func->utau(uground, fch);
        m_incflo->set_abl_friction_vels(uground, utau);
    }

    if (m_has_driving_dpdx) {
        const amrex::Real zh = m_abl_forcing->forcing_height();
        const amrex::Real vx = pa.line_velocity_xdir(zh);
        const amrex::Real vy = pa.line_velocity_ydir(zh);
        m_incflo->set_mean_abl_vel(vx, vy);
        m_abl_forcing->set_mean_velocities(vx, vy);
    }

    {
        // TODO: This should be handled by PlaneAveraging
        int output_interval = 1;
        amrex::ParmParse pp("amr");
        pp.query("line_plot_int", output_interval);

        if ((output_interval > 0) && (m_time.time_index() % output_interval == 0)) {
            pa.plot_line_text("line_plot.txt", m_time.time_index(), m_time.current_time());
        }
    }
}

} // namespace amr_wind

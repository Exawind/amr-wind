#include "amr-wind/physics/ConvectingTaylorVortex.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/utilities/trig_ops.H"

namespace amr_wind {
namespace ctv {

namespace {

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real UExact::operator()(
    const amrex::Real u0,
    const amrex::Real v0,
    const amrex::Real omega,
    const amrex::Real x,
    const amrex::Real y,
    const amrex::Real t) const
{
    return u0 - std::cos(utils::pi() * (x - u0 * t)) *
                    std::sin(utils::pi() * (y - v0 * t)) *
                    std::exp(-2.0 * omega * t);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real VExact::operator()(
    const amrex::Real u0,
    const amrex::Real v0,
    const amrex::Real omega,
    const amrex::Real x,
    const amrex::Real y,
    const amrex::Real t) const
{
    return v0 + std::sin(utils::pi() * (x - u0 * t)) *
                    std::cos(utils::pi() * (y - v0 * t)) *
                    std::exp(-2.0 * omega * t);
}

} // namespace

ConvectingTaylorVortex::ConvectingTaylorVortex(const CFDSim& sim)
    : m_time(sim.time())
    , m_repo(sim.repo())
    , m_mesh(sim.mesh())
    , m_velocity(sim.repo().get_field("velocity"))
    , m_density(sim.repo().get_field("density"))
{
    amrex::ParmParse pp("CTV");
    pp.query("density", m_rho);
    pp.query("u0", m_u0);
    pp.query("v0", m_v0);
    {
        amrex::Real nu;
        amrex::ParmParse pp("transport");
        pp.query("viscosity", nu);
        m_omega = utils::pi() * nu;
    }
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f;
        f.open(m_output_fname.c_str());
        f << std::setw(m_w) << "time" << std::setw(m_w) << "L2_u"
          << std::setw(m_w) << "L2_v" << std::endl;
        f.close();
    }
}

/** Initialize the velocity and density fields at the beginning of the
 *  simulation.
 */
void ConvectingTaylorVortex::initialize_fields(
    int level, const amrex::Geometry& geom)
{
    using namespace utils;

    const auto u0 = m_u0;
    const auto v0 = m_v0;
    const auto omega = m_omega;

    auto& velocity = m_velocity(level);
    auto& density = m_density(level);

    density.setVal(m_rho);

    UExact u_exact;
    VExact v_exact;

    for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();

        const auto& dx = geom.CellSizeArray();
        const auto& problo = geom.ProbLoArray();
        auto vel = velocity.array(mfi);

        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                vel(i, j, k, 0) = u_exact(u0, v0, omega, x, y, 0.0);
                vel(i, j, k, 1) = v_exact(u0, v0, omega, x, y, 0.0);
                vel(i, j, k, 2) = 0.0;
            });
    }
}

template <typename T>
amrex::Real ConvectingTaylorVortex::compute_error(const Field& field)
{

    amrex::Real error = 0.0;
    const amrex::Real time = m_time.current_time();
    const auto u0 = m_u0;
    const auto v0 = m_v0;
    const auto omega = m_omega;
    T f_exact;
    const auto comp = f_exact.m_comp;

    const int nlevels = m_repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {

        amrex::iMultiFab level_mask;
        if (lev < nlevels - 1) {
            level_mask = makeFineMask(
                m_mesh.boxArray(lev), m_mesh.DistributionMap(lev),
                m_mesh.boxArray(lev + 1), amrex::IntVect(2), 1, 0);
        } else {
            level_mask.define(
                m_mesh.boxArray(lev), m_mesh.DistributionMap(lev), 1, 0,
                amrex::MFInfo());
            level_mask.setVal(1);
        }

        const auto& dx = m_mesh.Geom(lev).CellSizeArray();
        const auto& problo = m_mesh.Geom(lev).ProbLoArray();
        const amrex::Real cell_vol = dx[0] * dx[1] * dx[2];

        const auto& fld = field(lev);
        error += amrex::ReduceSum(
            fld, level_mask, 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& fld_arr,
                amrex::Array4<int const> const& mask_arr) -> amrex::Real {
                amrex::Real err_fab = 0.0;

                amrex::Loop(bx, [=, &err_fab](int i, int j, int k) noexcept {
                    const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                    const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                    const amrex::Real u = fld_arr(i, j, k, comp);
                    const amrex::Real u_exact =
                        f_exact(u0, v0, omega, x, y, time);
                    err_fab += cell_vol * mask_arr(i, j, k) * (u - u_exact) *
                               (u - u_exact);
                });
                return err_fab;
            });
    }

    amrex::ParallelDescriptor::ReduceRealSum(error);

    const amrex::Real total_vol = m_mesh.Geom(0).ProbDomain().volume();
    return std::sqrt(error / total_vol);
}

void ConvectingTaylorVortex::output_error()
{
    const amrex::Real u_err = compute_error<UExact>(m_velocity);
    const amrex::Real v_err = compute_error<VExact>(m_velocity);

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f;
        f.open(m_output_fname.c_str(), std::ios_base::app);
        f << std::setprecision(12) << std::setw(m_w) << m_time.new_time()
          << std::setw(m_w) << u_err << std::setw(m_w) << v_err << std::endl;
        f.close();
    }
}

void ConvectingTaylorVortex::post_init_actions() { output_error(); }

void ConvectingTaylorVortex::post_advance_work() { output_error(); }

} // namespace ctv
} // namespace amr_wind

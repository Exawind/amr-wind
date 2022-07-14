#include "amr-wind/physics/ScalarAdvection.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real Gaussian::operator()(
    const amrex::Real x,
    const amrex::Real center,
    const amrex::Real amplitude,
    const amrex::Real width,
    const amrex::Real eta) const
{
    if (std::abs(x - center) < 12 * width/2) {
        return amplitude * std::exp( - std::pow(x - center, 2) / (2 * std::pow(width, 2)));
    } else {
        return 0.0;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real TopHat::operator()(
    const amrex::Real x,
    const amrex::Real center,
    const amrex::Real amplitude,
    const amrex::Real width,
    const amrex::Real eta) const
{
    if (std::abs(x - center) < width/2) {
        return amplitude;
    } else {
        return 0.0;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real Yalla2021::operator()(
    const amrex::Real x,
    const amrex::Real center,
    const amrex::Real amplitude,
    const amrex::Real width,
    const amrex::Real eta) const
{
    if (std::abs(x - center) < 12 * width/2) {
        return std::cos(eta * x) * amplitude * std::exp( - std::pow(x - center, 2) / (2 * std::pow(width, 2)));
    } else {
        return 0.0;
    }
}

ScalarAdvection::ScalarAdvection(CFDSim& sim)
    : m_sim(sim)
    , m_repo(sim.repo())
    , m_velocity(sim.repo().get_field("velocity"))
    , m_density(sim.repo().get_field("density"))
{
    // Register temperature equation
    auto& teqn = sim.pde_manager().register_transport_pde("Temperature");

    // Defer getting temperature field until PDE has been registered
    m_scalar = &(teqn.fields().field);

    amrex::ParmParse pp_scalar_advection("scalaradvection");
    pp_scalar_advection.query("center", m_center);
    pp_scalar_advection.query("amplitude", m_amplitude);
    pp_scalar_advection.query("width", m_width);
    pp_scalar_advection.query("eta", m_eta);
    pp_scalar_advection.query("shape", m_shape);

    amrex::ParmParse pp_incflo("incflo");
    pp_incflo.query("density", m_rho);
}

/** Initialize the velocity and temperature fields at the beginning of the
 *  simulation.
 */
void ScalarAdvection::initialize_fields(int level, const amrex::Geometry& geom)
{
    auto& velocity = m_velocity(level);
    auto& density = m_density(level);
    density.setVal(m_rho);

    for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        auto vel = velocity.array(mfi);

        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            vel(i, j, k, 0) = 1.0;
            vel(i, j, k, 1) = 0.0;
            vel(i, j, k, 2) = 0.0;
        });
    }
}

void ScalarAdvection::post_init_actions()
{
    if (m_sim.time().time_index() > 0) {
        return;
    }
    
    if (m_shape == "gaussian") {
        initialize_scalar(Gaussian());
    } else if (m_shape == "tophat") {
        initialize_scalar(TopHat());
    } else if (m_shape == "yalla2021") {
        initialize_scalar(Yalla2021());
    }
}

template <typename Shape>
void ScalarAdvection::initialize_scalar(const Shape& scalar_function)
{
    const amrex::Real center = m_center;
    const amrex::Real amplitude = m_amplitude;
    const amrex::Real width = m_width;
    const amrex::Real eta = m_eta;

    for (int level = 0; level <= m_repo.mesh().finestLevel(); ++level) {

        const auto& problo = m_repo.mesh().Geom(level).ProbLoArray();
        auto& scalar = (*m_scalar)(level);

        for (amrex::MFIter mfi(m_velocity(level)); mfi.isValid(); ++mfi) {
            const auto& dx = m_repo.mesh().Geom(level).CellSizeArray();
            const auto& nbx = mfi.nodaltilebox();
            auto sclr = scalar.array(mfi);

            amrex::ParallelFor(nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                sclr(i, j, k, 0) = scalar_function(x, center, amplitude, width, eta);
            });
        }
    }

}


} // namespace amr_wind

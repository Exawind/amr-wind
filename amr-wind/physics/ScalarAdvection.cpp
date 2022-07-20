#include "amr-wind/physics/ScalarAdvection.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"
#include "AMReX_MultiFabUtil.H"

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
    , m_time(sim.time())
    , m_repo(sim.repo())
    , m_velocity(sim.repo().get_field("velocity"))
    , m_density(sim.repo().get_field("density"))
{
    // Register temperature equation
    auto& teqn = sim.pde_manager().register_transport_pde("Temperature");

    // Defer getting temperature field until PDE has been registered
    m_scalar = &(teqn.fields().field);

    amrex::ParmParse pp_scalar_advection("scalaradvection");
    pp_scalar_advection.query("convective_speed", m_convective_speed);
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
    const amrex::Real convective_speed = m_convective_speed;
    auto& velocity = m_velocity(level);
    auto& density = m_density(level);
    density.setVal(m_rho);

    for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        auto vel = velocity.array(mfi);

        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            vel(i, j, k, 0) = convective_speed;
            vel(i, j, k, 1) = 0.0;
            vel(i, j, k, 2) = 0.0;
        });
    }
}

void ScalarAdvection::post_init_actions()
{
    if (m_time.time_index() > 0) {
        return;
    }
    
    if (m_shape == "gaussian") {
        initialize_scalar(Gaussian());
    } else if (m_shape == "tophat") {
        initialize_scalar(TopHat());
    } else if (m_shape == "yalla2021") {
        initialize_scalar(Yalla2021());
    }

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f;
        f.open(m_output_fname.c_str(), std::ofstream::out | std::ofstream::trunc);
        f.close();
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

        for (amrex::MFIter mfi(scalar); mfi.isValid(); ++mfi) {
            const auto& dx = m_repo.mesh().Geom(level).CellSizeArray();
            const auto& nbx = mfi.nodaltilebox();
            auto scalar_arr = scalar.array(mfi);

            amrex::ParallelFor(nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                scalar_arr(i, j, k, 0) = scalar_function(x, center, amplitude, width, eta);
            });
        }
    }

}

template <typename Shape>
amrex::Vector<amrex::Real> ScalarAdvection::compute_error(const Shape& scalar_function)
{
    const amrex::Real t = m_time.current_time();
    const amrex::Real convective_speed = m_convective_speed;
    const amrex::Real center = m_center;
    const amrex::Real amplitude = m_amplitude;
    const amrex::Real width = m_width;
    const amrex::Real eta = m_eta;
    const amrex::Real convected_center = center + convective_speed * t;
    const int nlevels = m_repo.num_active_levels();
    const amrex::Real total_vol = m_repo.mesh().Geom(0).ProbDomain().volume();

    amrex::Vector<amrex::Real> err(nlevels,0.0);

    for (int level = 0; level < nlevels; ++level) {
        amrex::Real err_lev = 0.0;
        amrex::iMultiFab level_mask;
        if (level < nlevels - 1) {
            level_mask = makeFineMask(
                m_repo.mesh().boxArray(level), m_repo.mesh().DistributionMap(level),
                m_repo.mesh().boxArray(level + 1), amrex::IntVect(2), 1, 0);
        } else {
            level_mask.define(
                m_repo.mesh().boxArray(level), m_repo.mesh().DistributionMap(level), 1, 0,
                amrex::MFInfo());
            level_mask.setVal(1);
        }

        const auto& dx = m_repo.mesh().Geom(level).CellSizeArray();
        const auto& problo = m_repo.mesh().Geom(level).ProbLoArray();
        const amrex::Real cell_vol = dx[0] * dx[1] * dx[2];

        const auto& scalar = (*m_scalar)(level);
        for (amrex::MFIter mfi(scalar); mfi.isValid(); ++mfi) {
            const auto& vbx = mfi.validbox();
            const auto& scalar_arr = scalar.array(mfi);
            const auto& mask_arr = level_mask.array(mfi);

            amrex::Real err_fab = 0.0;
            amrex::LoopOnCpu(vbx, [=, &err_fab](int i, int j, int k) noexcept {
                const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                const amrex::Real s = scalar_arr(i, j, k, 0);
                const amrex::Real s_exact = scalar_function(x, convected_center, amplitude, width, eta);
                err_fab += cell_vol * mask_arr(i, j, k) * (s - s_exact) *
                           (s - s_exact);
            });
            err_lev += err_fab;
        }
        amrex::ParallelDescriptor::ReduceRealSum(err_lev);
        err[level] = std::sqrt(err_lev / total_vol);
    }

    return err;
}

void ScalarAdvection::post_advance_work()
{
    const int nlevels = m_repo.num_active_levels();
    amrex::Vector<amrex::Real> err(nlevels,0.0);

    if (m_shape == "gaussian") {
        err = compute_error(Gaussian());
    } else if (m_shape == "tophat") {
        err = compute_error(TopHat());
    } else if (m_shape == "yalla2021") {
        err = compute_error(Yalla2021());
    }

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f;
        f.open(m_output_fname.c_str(), std::ios_base::app);
        f << std::setprecision(12) << std::setw(m_w) << m_time.current_time()
        << std::setw(m_w);
        for(auto i: err){
            f << i << std::setw(m_w);
        }
        f << std::endl;
        f.close();
    }
}

} // namespace amr_wind

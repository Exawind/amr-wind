#include "amr-wind/physics/ScalarAdvection.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"
#include "AMReX_MultiFabUtil.H"

namespace amr_wind {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
GaussianPulseFV::operator()(
    const amrex::Real x,
    const amrex::Real /*unused*/,
    const amrex::Real dx,
    const amrex::Real /*unused*/,
    const amrex::Real x0,
    const amrex::Real /*unused*/,
    const amrex::Real amplitude,
    const amrex::Real x_width,
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/) const
{
    amrex::Real val = 0.0;
    if (std::abs(x - x0) < 6 * x_width) {
        val = sqrt(utils::pi() / 2) * amplitude * x_width *
              (std::erf((x - x0 + dx / 2) / (sqrt(2) * x_width)) -
               std::erf((x - x0 - dx / 2) / (sqrt(2) * x_width))) /
              dx;
    }
    return val;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
TwoDimGaussianPulseFV::operator()(
    const amrex::Real x,
    const amrex::Real y,
    const amrex::Real dx,
    const amrex::Real dy,
    const amrex::Real x0,
    const amrex::Real y0,
    const amrex::Real amplitude,
    const amrex::Real x_width,
    const amrex::Real y_width,
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/) const
{
    amrex::Real val = 0.0;
    if (std::abs(x - x0) < 6 * x_width && std::abs(y - y0) < 6 * y_width) {
        val = utils::pi() / 2 * amplitude * x_width * y_width *
              (std::erf((x - x0 + dx / 2) / (sqrt(2) * x_width)) -
               std::erf((x - x0 - dx / 2) / (sqrt(2) * x_width))) *
              (std::erf((y - y0 + dy / 2) / (sqrt(2) * y_width)) -
               std::erf((y - y0 - dy / 2) / (sqrt(2) * y_width))) /
              dx / dy;
    }
    return val;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real SquarePulseFV::operator()(
    const amrex::Real x,
    const amrex::Real /*unused*/,
    const amrex::Real dx,
    const amrex::Real /*unused*/,
    const amrex::Real x0,
    const amrex::Real /*unused*/,
    const amrex::Real amplitude,
    const amrex::Real x_width,
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/) const
{
    amrex::Real val = 0.0;
    if (std::abs(std::abs(x - x0) - x_width / 2) < dx / 2) {
        val = amplitude * (x_width / 2 - std::abs(x - x0) + dx / 2) / dx;
    } else if (std::abs(x - x0) < x_width / 2) {
        val = amplitude;
    }
    return val;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
GaussianWavePacketFV::operator()(
    const amrex::Real x,
    const amrex::Real /*unused*/,
    const amrex::Real dx,
    const amrex::Real /*unused*/,
    const amrex::Real x0,
    const amrex::Real /*unused*/,
    const amrex::Real amplitude,
    const amrex::Real x_width,
    const amrex::Real /*unused*/,
    const amrex::Real x_wavenumber,
    const amrex::Real /*unused*/) const
{
    const GaussianWavePacket pointwise_function;
    // Gauss-Legendre quadrature
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx_i = {
        AMREX_D_DECL(-0.7745966692, 0, 0.7745966692)};
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> w_i = {
        AMREX_D_DECL(0.5555555556, 0.8888888889, 0.5555555556)};
    amrex::Real cell_integral = 0.0;
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        cell_integral =
            cell_integral + w_i[i] * pointwise_function(
                                         x + dx_i[i] * 0.5 * dx, x0, amplitude,
                                         x_width, x_wavenumber);
    }
    return cell_integral / 2;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
GaussianWavePacket::operator()(
    const amrex::Real x,
    const amrex::Real x0,
    const amrex::Real amplitude,
    const amrex::Real x_width,
    const amrex::Real x_wavenumber) const
{
    amrex::Real val = 0.0;
    if (std::abs(x - x0) < 6 * x_width) {
        val = amplitude * std::cos(x_wavenumber * x) *
              std::exp(-std::pow(x - x0, 2) / (2 * std::pow(x_width, 2)));
    }
    return val;
}

ScalarAdvection::ScalarAdvection(CFDSim& sim)
    : m_time(sim.time())
    , m_repo(sim.repo())
    , m_velocity(sim.repo().get_field("velocity"))
    , m_density(sim.repo().get_field("density"))
{
    // Register temperature equation
    auto& teqn = sim.pde_manager().register_transport_pde("Temperature");

    // Defer getting temperature field until PDE has been registered
    m_scalar = &(teqn.fields().field);

    amrex::ParmParse pp_scalar_advection("scalaradvection");
    pp_scalar_advection.query("u", m_u);
    pp_scalar_advection.query("v", m_v);
    pp_scalar_advection.query("x0", m_x0);
    pp_scalar_advection.query("y0", m_y0);
    pp_scalar_advection.query("amplitude", m_amplitude);
    pp_scalar_advection.query("x_width", m_x_width);
    pp_scalar_advection.query("y_width", m_y_width);
    pp_scalar_advection.query("x_wavenumber", m_x_wavenumber);
    pp_scalar_advection.query("y_wavenumber", m_y_wavenumber);
    pp_scalar_advection.query("shape", m_shape);
    pp_scalar_advection.query("output_fname", m_output_fname);

    amrex::ParmParse pp_incflo("incflo");
    pp_incflo.query("density", m_rho);
}

/** Initialize the velocity and temperature fields at the beginning of the
 *  simulation.
 */
void ScalarAdvection::initialize_fields(
    int level, const amrex::Geometry& /*geom*/)
{
    const amrex::Real u = m_u;
    const amrex::Real v = m_v;
    auto& velocity = m_velocity(level);
    auto& density = m_density(level);
    density.setVal(m_rho);

    for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        auto vel = velocity.array(mfi);

        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                vel(i, j, k, 0) = u;
                vel(i, j, k, 1) = v;
                vel(i, j, k, 2) = 0.0;
            });
    }
}

void ScalarAdvection::post_init_actions()
{
    if (m_time.time_index() > 0) {
        return;
    }

    if (m_shape == "gaussianpulse") {
        initialize_scalar(GaussianPulseFV());
    } else if (m_shape == "squarepulse") {
        initialize_scalar(SquarePulseFV());
    } else if (m_shape == "gaussianwavepacket") {
        initialize_scalar(GaussianWavePacketFV());
    } else if (m_shape == "twodimgaussianpulse") {
        initialize_scalar(TwoDimGaussianPulseFV());
    }

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f;
        f.open(
            m_output_fname.c_str(), std::ofstream::out | std::ofstream::trunc);
        f.close();
    }
}

template <typename Shape>
void ScalarAdvection::initialize_scalar(const Shape& scalar_function)
{
    const amrex::Real x0 = m_x0;
    const amrex::Real y0 = m_y0;
    const amrex::Real amplitude = m_amplitude;
    const amrex::Real x_width = m_x_width;
    const amrex::Real y_width = m_y_width;
    const amrex::Real x_wavenumber = m_x_wavenumber;
    const amrex::Real y_wavenumber = m_y_wavenumber;

    for (int level = 0; level <= m_repo.mesh().finestLevel(); ++level) {

        const auto& problo = m_repo.mesh().Geom(level).ProbLoArray();
        auto& scalar = (*m_scalar)(level);

        for (amrex::MFIter mfi(scalar); mfi.isValid(); ++mfi) {
            const auto& dx = m_repo.mesh().Geom(level).CellSizeArray();
            const auto& nbx = mfi.nodaltilebox();
            auto scalar_arr = scalar.array(mfi);

            amrex::ParallelFor(
                nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                    const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                    scalar_arr(i, j, k, 0) = scalar_function(
                        x, y, dx[0], dx[1], x0, y0, amplitude, x_width, y_width,
                        x_wavenumber, y_wavenumber);
                });
        }
    }
}

template <typename Shape>
amrex::Vector<amrex::Real>
ScalarAdvection::compute_error(const Shape& scalar_function)
{
    const amrex::Real t = m_time.new_time();
    const amrex::Real x_conv_dist = m_u * t;
    const amrex::Real y_conv_dist = m_v * t;
    const amrex::Real x0 = m_x0;
    const amrex::Real y0 = m_y0;
    const amrex::Real amplitude = m_amplitude;
    const amrex::Real x_width = m_x_width;
    const amrex::Real y_width = m_y_width;
    const amrex::Real x_wavenumber = m_x_wavenumber;
    const amrex::Real y_wavenumber = m_y_wavenumber;
    const int nlevels = m_repo.num_active_levels();
    const amrex::Real total_vol = m_repo.mesh().Geom(0).ProbDomain().volume();

    amrex::Vector<amrex::Real> err(nlevels + 1, 0.0);

    for (int level = 0; level < nlevels; ++level) {
        amrex::Real err_lev = 0.0;
        amrex::iMultiFab level_mask;
        if (level < nlevels - 1) {
            level_mask = makeFineMask(
                m_repo.mesh().boxArray(level),
                m_repo.mesh().DistributionMap(level),
                m_repo.mesh().boxArray(level + 1), amrex::IntVect(2), 1, 0);
        } else {
            level_mask.define(
                m_repo.mesh().boxArray(level),
                m_repo.mesh().DistributionMap(level), 1, 0, amrex::MFInfo());
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
                const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                const amrex::Real s = scalar_arr(i, j, k, 0);
                const amrex::Real s_exact = scalar_function(
                    x - x_conv_dist, y - y_conv_dist, dx[0], dx[1], x0, y0,
                    amplitude, x_width, y_width, x_wavenumber, y_wavenumber);
                err_fab += cell_vol * mask_arr(i, j, k) * (s - s_exact) *
                           (s - s_exact);
            });
            err_lev += err_fab;
        }
        amrex::ParallelDescriptor::ReduceRealSum(err_lev);
        err[level] = err_lev;
        err[nlevels] += err_lev;
    }

    for (int level = 0; level < nlevels + 1; ++level) {
        err[level] = std::sqrt(err[level] / total_vol);
    }

    return err;
}

void ScalarAdvection::post_advance_work()
{
    const int nlevels = m_repo.num_active_levels();
    amrex::Vector<amrex::Real> err(nlevels + 1, 0.0);

    if (m_shape == "gaussianpulse") {
        err = compute_error(GaussianPulseFV());
    } else if (m_shape == "squarepulse") {
        err = compute_error(SquarePulseFV());
    } else if (m_shape == "gaussianwavepacket") {
        err = compute_error(GaussianWavePacketFV());
    } else if (m_shape == "twodimgaussianpulse") {
        err = compute_error(TwoDimGaussianPulseFV());
    }

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f;
        f.open(m_output_fname.c_str(), std::ios_base::app);
        f << std::setprecision(12) << std::setw(m_w) << m_time.new_time()
          << std::setw(m_w);
        for (auto i : err) {
            f << i << std::setw(m_w);
        }
        f << std::endl;
        f.close();
    }
}

} // namespace amr_wind

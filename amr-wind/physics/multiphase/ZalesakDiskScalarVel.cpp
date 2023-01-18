#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/physics/multiphase/ZalesakDiskScalarVel.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"
#include "AMReX_MultiFabUtil.H"
#include "amr-wind/fvm/gradient.H"
#include "amr-wind/core/field_ops.H"

namespace amr_wind::zds {

namespace {
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real SCexact::operator()(
    amrex::Real xc0,
    amrex::Real yc0,
    const amrex::Real zc0,
    const amrex::Real x,
    const amrex::Real y,
    const amrex::Real z,
    const amrex::Real radius,
    const amrex::Real alpha) const
{
    // Adjust initial x and y to be relative to center of domain
    // (assuming 0,1 boundary)
    xc0 = xc0 - 0.5;
    yc0 = yc0 - 0.5;
    // Calculate new center based on time
    amrex::Real xc = 0.5 + std::cos(alpha) * xc0 - std::sin(alpha) * yc0;
    amrex::Real yc = 0.5 + std::sin(alpha) * xc0 + std::cos(alpha) * yc0;
    amrex::Real zc = zc0;

    // Normalized distance from current center
    amrex::Real dnorm =
        std::sqrt(
            (x - xc) * (x - xc) + (y - yc) * (y - yc) + (z - zc) * (z - zc)) /
        radius;
    // Scalar distribution
    return amrex::min(1.0, amrex::max(0.0, 1.5 - dnorm));
}
} // namespace

ZalesakDiskScalarVel::ZalesakDiskScalarVel(CFDSim& sim)
    : m_sim(sim)
    , m_velocity(sim.repo().get_field("velocity"))
    , m_levelset(sim.repo().get_field("levelset"))
    , m_density(sim.repo().get_field("density"))
{
    amrex::ParmParse pp(identifier());
    pp.queryarr("location", m_loc, 0, AMREX_SPACEDIM);
    pp.query("radius", m_radius);
    pp.query("period", m_TT);
    pp.query("error_log_file", m_output_fname);
    amrex::ParmParse pinc("incflo");
    pinc.add("prescribe_velocity", true);
}

/** Initialize the velocity and levelset fields at the beginning of the
 *  simulation.
 *
 *  \sa amr_wind::ZalesakDiskScalarVelFieldInit
 */
void ZalesakDiskScalarVel::initialize_fields(
    int level, const amrex::Geometry& geom)
{
    auto& velocity = m_velocity(level);
    auto& levelset = m_levelset(level);
    auto& density = m_density(level);
    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();

    const auto& mphase = m_sim.physics_manager().get<MultiPhase>();
    const amrex::Real rho1 = mphase.rho1();
    const amrex::Real rho2 = mphase.rho2();

    auto& u_mac = m_sim.repo().get_field("u_mac")(level);
    auto& v_mac = m_sim.repo().get_field("v_mac")(level);
    auto& w_mac = m_sim.repo().get_field("w_mac")(level);

    const amrex::Real xc = m_loc[0];
    const amrex::Real yc = m_loc[1];
    const amrex::Real zc = m_loc[2];
    const amrex::Real radius = m_radius;
    const amrex::Real TT = m_TT;
    const amrex::Real width = m_width;
    const amrex::Real depth = m_depth;

    for (amrex::MFIter mfi(levelset); mfi.isValid(); ++mfi) {
        const auto& gbx = mfi.growntilebox(1);
        auto uf = u_mac.array(mfi);
        auto vf = v_mac.array(mfi);
        auto wf = w_mac.array(mfi);
        auto vel = velocity.array(mfi);
        auto phi = levelset.array(mfi);
        auto rho = density.array(mfi);
        amrex::ParallelFor(
            gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

                uf(i, j, k) = 2.0 * M_PI / TT * (0.5 - y);
                vf(i, j, k) = 2.0 * M_PI / TT * (x - 0.5);
                wf(i, j, k) = 0.0;

                vel(i, j, k, 1) = 0.0;
                vel(i, j, k, 2) = 0.0;
                // First define the sphere

                phi(i, j, k) =
                    radius - std::sqrt(
                                 (x - xc) * (x - xc) + (y - yc) * (y - yc) +
                                 (z - zc) * (z - zc));
                // then the slot
                amrex::Real eps = std::cbrt(dx[0] * dx[1] * dx[2]);
                if (y - yc <= radius && y - yc >= radius - depth &&
                    std::abs(x - xc) <= width &&
                    std::sqrt(
                        (x - xc) * (x - xc) + (y - yc) * (y - yc) +
                        (z - zc) * (z - zc)) < radius + eps) {
                    amrex::Real d1;
                    if (x > xc) {
                        d1 = std::abs(xc + width - x);
                    } else {
                        d1 = std::abs(xc - width - x);
                    }
                    const amrex::Real d2 = std::abs(y - (yc + radius - depth));
                    const amrex::Real min_dist = amrex::min(d1, d2);

                    phi(i, j, k) = -min_dist;
                }

                amrex::Real smooth_heaviside;
                eps = std::cbrt(2. * dx[0] * dx[1] * dx[2]);
                if (phi(i, j, k) > eps) {
                    smooth_heaviside = 1.0;
                } else if (phi(i, j, k) < -eps) {
                    smooth_heaviside = 0.;
                } else {
                    smooth_heaviside =
                        0.5 *
                        (1.0 + phi(i, j, k) / eps +
                         1.0 / M_PI * std::sin(phi(i, j, k) * M_PI / eps));
                }
                rho(i, j, k) =
                    rho1 * smooth_heaviside + rho2 * (1.0 - smooth_heaviside);

                // Set up scalar field with velocity in x
                amrex::Real dnorm =
                    std::sqrt(
                        (x - xc) * (x - xc) + (y - yc) * (y - yc) +
                        (z - zc) * (z - zc)) /
                    radius;
                // Set up scalar field with u velocity
                vel(i, j, k, 0) = amrex::min(1.0, amrex::max(0.0, 1.5 - dnorm));
            });
    }
    m_levelset.fillpatch(0.0);
    m_velocity.fillpatch(0.0);
    m_density.fillpatch(0.0);
    output_error();
}

void ZalesakDiskScalarVel::pre_advance_work()
{

    const int nlevels = m_sim.repo().num_active_levels();
    const auto& geom = m_sim.mesh().Geom();

    // Overriding the velocity field
    for (int lev = 0; lev < nlevels; ++lev) {
        auto& u_mac = m_sim.repo().get_field("u_mac")(lev);
        auto& v_mac = m_sim.repo().get_field("v_mac")(lev);
        auto& w_mac = m_sim.repo().get_field("w_mac")(lev);
        for (amrex::MFIter mfi(m_velocity(lev)); mfi.isValid(); ++mfi) {
            const auto& gbx = mfi.growntilebox(1);
            const auto& dx = geom[lev].CellSizeArray();
            const auto& problo = geom[lev].ProbLoArray();
            const amrex::Real TT = m_TT;
            auto uf = u_mac.array(mfi);
            auto vf = v_mac.array(mfi);
            auto wf = w_mac.array(mfi);
            amrex::ParallelFor(
                gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                    const amrex::Real y = problo[1] + (j + 0.5) * dx[1];

                    uf(i, j, k) = 2.0 * M_PI / TT * (0.5 - y);
                    vf(i, j, k) = 2.0 * M_PI / TT * (x - 0.5);
                    wf(i, j, k) = 0.0;
                });
        }
        u_mac.FillBoundary(geom[lev].periodicity());
        v_mac.FillBoundary(geom[lev].periodicity());
        w_mac.FillBoundary(geom[lev].periodicity());
    }
}

void ZalesakDiskScalarVel::post_advance_work() { output_error(); }

template <typename T>
amrex::Real ZalesakDiskScalarVel::compute_error(const Field& field)
{
    amrex::Real error = 0.0;
    const amrex::Real xc = m_loc[0];
    const amrex::Real yc = m_loc[1];
    const amrex::Real zc = m_loc[2];
    const amrex::Real radius = m_radius;
    const amrex::Real TT = m_TT;
    const amrex::Real time = m_sim.time().new_time();
    T f_exact;
    const auto comp = f_exact.m_comp;
    const auto mesh_mapping = m_sim.has_mesh_mapping();

    Field const* nu_coord_cc =
        mesh_mapping ? &(m_sim.repo().get_field("non_uniform_coord_cc"))
                     : nullptr;
    Field const* mesh_fac_cc =
        mesh_mapping
            ? &(m_sim.repo().get_mesh_mapping_field(amr_wind::FieldLoc::CELL))
            : nullptr;

    const int nlevels = m_sim.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {

        amrex::iMultiFab level_mask;
        if (lev < nlevels - 1) {
            level_mask = makeFineMask(
                m_sim.mesh().boxArray(lev), m_sim.mesh().DistributionMap(lev),
                m_sim.mesh().boxArray(lev + 1), amrex::IntVect(2), 1, 0);
        } else {
            level_mask.define(
                m_sim.mesh().boxArray(lev), m_sim.mesh().DistributionMap(lev),
                1, 0, amrex::MFInfo());
            level_mask.setVal(1);
        }

        if (m_sim.has_overset()) {
            for (amrex::MFIter mfi(field(lev)); mfi.isValid(); ++mfi) {
                const auto& vbx = mfi.validbox();

                const auto& iblank_arr =
                    m_sim.repo().get_int_field("iblank_cell")(lev).array(mfi);
                const auto& imask_arr = level_mask.array(mfi);
                amrex::ParallelFor(
                    vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        if (amrex::Math::abs(iblank_arr(i, j, k)) < 1) {
                            imask_arr(i, j, k) = 0;
                        }
                    });
            }
        }

        const auto& dx = m_sim.mesh().Geom(lev).CellSizeArray();
        const auto& prob_lo = m_sim.mesh().Geom(lev).ProbLoArray();

        const auto& fld = field(lev);
        auto const& fld_arr = fld.const_arrays();
        auto const& mask_arr = level_mask.const_arrays();
        amrex::MultiArray4<amrex::Real const> fac_arr =
            mesh_mapping ? ((*mesh_fac_cc)(lev).const_arrays())
                         : amrex::MultiArray4<amrex::Real const>();
        amrex::MultiArray4<amrex::Real const> nu_cc =
            mesh_mapping ? ((*nu_coord_cc)(lev).const_arrays())
                         : amrex::MultiArray4<amrex::Real const>();

        error += amrex::ParReduce(
            amrex::TypeList<amrex::ReduceOpSum>{},
            amrex::TypeList<amrex::Real>{}, fld, amrex::IntVect(0),
            [=] AMREX_GPU_HOST_DEVICE(int box_no, int i, int j, int k)
                -> amrex::GpuTuple<amrex::Real> {
                auto const& fld_bx = fld_arr[box_no];
                auto const& mask_bx = mask_arr[box_no];

                amrex::Real x = mesh_mapping ? (nu_cc[box_no](i, j, k, 0))
                                             : (prob_lo[0] + (i + 0.5) * dx[0]);
                amrex::Real y = mesh_mapping ? (nu_cc[box_no](i, j, k, 1))
                                             : (prob_lo[1] + (j + 0.5) * dx[1]);
                amrex::Real z = mesh_mapping ? (nu_cc[box_no](i, j, k, 2))
                                             : (prob_lo[2] + (k + 0.5) * dx[2]);
                amrex::Real fac_x =
                    mesh_mapping ? (fac_arr[box_no](i, j, k, 0)) : 1.0;
                amrex::Real fac_y =
                    mesh_mapping ? (fac_arr[box_no](i, j, k, 1)) : 1.0;
                amrex::Real fac_z =
                    mesh_mapping ? (fac_arr[box_no](i, j, k, 2)) : 1.0;

                const amrex::Real u = fld_bx(i, j, k, comp);
                const amrex::Real u_exact = f_exact(
                    xc, yc, zc, x, y, z, radius, time * 2.0 * M_PI / TT);
                const amrex::Real cell_vol =
                    dx[0] * fac_x * dx[1] * fac_y * dx[2] * fac_z;

                return cell_vol * mask_bx(i, j, k) * (u - u_exact) *
                       (u - u_exact);
            });
    }

    amrex::ParallelDescriptor::ReduceRealSum(error);

    const amrex::Real total_vol = m_sim.mesh().Geom(0).ProbDomain().volume();
    return std::sqrt(error / total_vol);
}

void ZalesakDiskScalarVel::output_error()
{
    // TODO: gradp analytical solution has not been adjusted for mesh mapping
    const amrex::Real SC_err = compute_error<zds::SCexact>(m_velocity);

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f;
        f.open(m_output_fname.c_str(), std::ios_base::app);
        f << std::setprecision(12) << std::setw(m_w) << m_sim.time().new_time()
          << std::setw(m_w) << SC_err << std::endl;
        f.close();
    }
}
} // namespace amr_wind::zds

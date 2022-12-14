#include "amr-wind/physics/ChannelFlow.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/utilities/DirectionSelector.H"

namespace amr_wind::channel_flow {

ChannelFlow::ChannelFlow(CFDSim& sim)
    : m_time(sim.time())
    , m_repo(sim.repo())
    , m_mesh(sim.mesh())
    , m_mesh_mapping(sim.has_mesh_mapping())
{
    {
        std::string turbulence_model;
        amrex::ParmParse pp("turbulence");
        pp.query("model", turbulence_model);
        if (turbulence_model == "Laminar") {
            m_laminar = true;
        }
    }
    {
        amrex::ParmParse pp("ChannelFlow");
        pp.query("flow_direction", m_mean_vel_dir);
        pp.query("normal_direction", m_norm_dir);
        pp.query("error_log_file", m_output_fname);

        if (m_laminar) {
            pp.query("density", m_rho);
            pp.query("Mean_Velocity", m_mean_vel);
        } else {
            pp.query("density", m_rho);
            pp.query("re_tau", m_re_tau);
            pp.query("tke0", m_tke0);
            pp.query("sdr0", m_sdr0);
        }
    }
    {
        amrex::ParmParse pp("transport");
        pp.query("viscosity", m_mu);
        // Assumes a boundary layer height of 1.0
        m_utau = m_mu * m_re_tau / (m_rho * 1.0);
        m_ytau = m_mu / (m_utau * m_rho);
    }
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f;
        f.open(m_output_fname.c_str());
        f << std::setw(m_w) << "time" << std::setw(m_w) << "L2_u" << std::endl;
        f.close();
    }
}

/** Initialize the velocity, density, tke and sdr fields at the beginning of the
 *  simulation.
 */
void ChannelFlow::initialize_fields(int level, const amrex::Geometry& geom)
{
    switch (m_norm_dir) {
    case 0:
        initialize_fields(level, geom, XDir(), 0);
        break;
    case 1:
        initialize_fields(level, geom, YDir(), 1);
        break;
    case 2:
        initialize_fields(level, geom, ZDir(), 2);
        break;
    default:
        amrex::Abort("axis must be equal to 1 or 2");
        break;
    }
}

template <typename IndexSelector>
void ChannelFlow::initialize_fields(
    int level,
    const amrex::Geometry& geom,
    const IndexSelector& idxOp,
    const int n_idx)
{
    const amrex::Real kappa = m_kappa;
    const amrex::Real y_tau = m_ytau;
    const amrex::Real utau = m_utau;
    auto& velocity_field = m_repo.get_field("velocity");
    auto& velocity = velocity_field(level);
    auto& density = m_repo.get_field("density")(level);
    density.setVal(m_rho);

    if (!m_laminar) {
        auto& tke = m_repo.get_field("tke")(level);
        auto& sdr = m_repo.get_field("sdr")(level);
        auto& walldist = m_repo.get_field("wall_dist")(level);
        tke.setVal(m_tke0);
        sdr.setVal(m_sdr0);

        for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
            const auto& vbx = mfi.validbox();

            const auto& dx = geom.CellSizeArray();
            const auto& problo = geom.ProbLoArray();
            auto vel = velocity.array(mfi);
            auto wd = walldist.array(mfi);

            amrex::ParallelFor(
                vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const int n_ind = idxOp(i, j, k);
                    amrex::Real h = problo[n_idx] + (n_ind + 0.5) * dx[n_idx];
                    if (h > 1.0) {
                        h = 2.0 - h;
                    }
                    wd(i, j, k) = h;
                    const amrex::Real hp = h / y_tau;
                    vel(i, j, k, 0) =
                        utau * (1. / kappa * std::log1p(kappa * hp) +
                                7.8 * (1.0 - std::exp(-hp / 11.0) -
                                       (hp / 11.0) * std::exp(-hp / 3.0)));
                    vel(i, j, k, 1) = 0.0;
                    vel(i, j, k, 2) = 0.0;
                });
        }
    } else {
        velocity.setVal(0.0);
        velocity.setVal(
            m_mean_vel, m_mean_vel_dir, 1,
            velocity_field.num_grow()[m_mean_vel_dir]);
    }
}

template <typename IndexSelector>
amrex::Real ChannelFlow::compute_error(const IndexSelector& idxOp)
{
    amrex::Real error = 0.0;
    const auto flow_dir = m_mean_vel_dir;
    const auto norm_dir = m_norm_dir;
    const auto mu = m_mu;
    const auto mesh_mapping = m_mesh_mapping;

    amrex::Real dpdx = 0;
    {
        amrex::Vector<amrex::Real> body_force{{0.0, 0.0, 0.0}};
        amrex::ParmParse pp("BodyForce");
        pp.queryarr("magnitude", body_force, 0, AMREX_SPACEDIM);
        dpdx = body_force[flow_dir];
    }

    const auto& problo = m_mesh.Geom(0).ProbLoArray();
    amrex::Real ht = 0.0;
    {
        amrex::Vector<amrex::Real> probhi_physical{{0.0, 0.0, 0.0}};
        amrex::ParmParse pp("geometry");
        if (pp.contains("prob_hi_physical")) {
            pp.getarr("prob_hi_physical", probhi_physical);
        } else {
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                probhi_physical[d] = m_mesh.Geom(0).ProbHiArray()[d];
            }
        }
        ht = probhi_physical[norm_dir] - problo[norm_dir];
    }

    const auto& velocity = m_repo.get_field("velocity");
    Field const* nu_coord_cc =
        mesh_mapping ? &(m_repo.get_field("non_uniform_coord_cc")) : nullptr;
    Field const* mesh_fac_cc =
        mesh_mapping
            ? &(m_repo.get_mesh_mapping_field(amr_wind::FieldLoc::CELL))
            : nullptr;

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
        const auto& prob_lo = m_mesh.Geom(lev).ProbLoArray();

        const auto& vel = velocity(lev);
        auto const& vel_arr = vel.const_arrays();
        auto const& mask_arr = level_mask.const_arrays();
        amrex::MultiArray4<amrex::Real const> fac_arr =
            mesh_mapping ? ((*mesh_fac_cc)(lev).const_arrays())
                         : amrex::MultiArray4<amrex::Real const>();
        amrex::MultiArray4<amrex::Real const> nu_cc =
            mesh_mapping ? ((*nu_coord_cc)(lev).const_arrays())
                         : amrex::MultiArray4<amrex::Real const>();

        error += amrex::ParReduce(
            amrex::TypeList<amrex::ReduceOpSum>{},
            amrex::TypeList<amrex::Real>{}, vel, amrex::IntVect(0),
            [=] AMREX_GPU_HOST_DEVICE(int box_no, int i, int j, int k)
                -> amrex::GpuTuple<amrex::Real> {
                auto const& vel_bx = vel_arr[box_no];
                auto const& mask_bx = mask_arr[box_no];

                amrex::Real y = mesh_mapping
                                    ? (nu_cc[box_no](i, j, k, norm_dir))
                                    : (prob_lo[norm_dir] +
                                       (idxOp(i, j, k) + 0.5) * dx[norm_dir]);
                amrex::Real fac_x =
                    mesh_mapping ? (fac_arr[box_no](i, j, k, 0)) : 1.0;
                amrex::Real fac_y =
                    mesh_mapping ? (fac_arr[box_no](i, j, k, 1)) : 1.0;
                amrex::Real fac_z =
                    mesh_mapping ? (fac_arr[box_no](i, j, k, 2)) : 1.0;

                const amrex::Real u = vel_bx(i, j, k, flow_dir);
                const amrex::Real u_exact =
                    1 / (2 * mu) * -dpdx * (y * y - y * ht);

                const amrex::Real cell_vol =
                    dx[0] * fac_x * dx[1] * fac_y * dx[2] * fac_z;

                return cell_vol * mask_bx(i, j, k) * (u - u_exact) *
                       (u - u_exact);
            });
    }

    amrex::ParallelDescriptor::ReduceRealSum(error);

    const amrex::Real total_vol = m_mesh.Geom(0).ProbDomain().volume();
    return std::sqrt(error / total_vol);
}

void ChannelFlow::output_error()
{
    // analytical solution exists only for flow in streamwise direction
    amrex::Real u_err = 0.0;
    switch (m_norm_dir) {
    case 0:
        u_err = compute_error(XDir());
        break;
    case 1:
        u_err = compute_error(YDir());
        break;
    case 2:
        u_err = compute_error(ZDir());
        break;
    default:
        amrex::Abort("axis must be equal to 1 or 2");
        break;
    }

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f;
        f.open(m_output_fname.c_str(), std::ios_base::app);
        f << std::setprecision(12) << std::setw(m_w) << m_time.new_time()
          << std::setw(m_w) << u_err << std::endl;
        f.close();
    }
}

void ChannelFlow::post_init_actions()
{
    if (m_laminar) {
        output_error();
    }
}

void ChannelFlow::post_advance_work()
{
    if (m_laminar) {
        output_error();
    }
}

} // namespace amr_wind::channel_flow

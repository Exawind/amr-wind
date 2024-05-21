#include "amr-wind/physics/BurggrafFlow.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParmParse.H"
#include "AMReX_ParReduce.H"
#include "amr-wind/utilities/trig_ops.H"

namespace amr_wind::burggraf {

namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
UExact::operator()(const amrex::Real x, const amrex::Real y) const
{
    return 8 * (std::pow(x, 4) - 2 * std::pow(x, 3) + std::pow(x, 2)) *
           (4 * std::pow(y, 3) - 2 * y);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
VExact::operator()(const amrex::Real x, const amrex::Real y) const
{
    return -8 * (4 * std::pow(x, 3) - 6 * std::pow(x, 2) + 2 * x) *
           (std::pow(y, 4) - std::pow(y, 2));
}

} // namespace

BurggrafFlow::BurggrafFlow(const CFDSim& sim)
    : m_time(sim.time())
    , m_repo(sim.repo())
    , m_mesh(sim.mesh())
    , m_velocity(sim.repo().get_field("velocity"))
    , m_density(sim.repo().get_field("density"))
    , m_source(sim.repo().declare_field("bf_src_term", 3))
{
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f;
        f.open(m_output_fname.c_str());
        f << std::setw(m_w) << "time" << std::setw(m_w) << "L2_u"
          << std::setw(m_w) << "L2_v" << std::setw(m_w) << "L2_w" << std::endl;
        f.close();
    }
}

void BurggrafFlow::initialize_fields(int level, const amrex::Geometry& geom)
{
    using namespace utils;

    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();

    auto& velocity = m_velocity(level);
    auto& density = m_density(level);
    auto& source = m_source(level);

    amrex::Real nu;
    amrex::ParmParse pp("transport");
    pp.query("viscosity", nu);
    amrex::Real Re = 1 / nu;

    density.setVal(m_rho);

    UExact u_exact;
    VExact v_exact;
    for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();

        auto vel = velocity.array(mfi);
        auto src = source.array(mfi);

        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // compute the source term
                const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
                const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];

                vel(i, j, k, 0) = u_exact(x, y);
                vel(i, j, k, 1) = v_exact(x, y);
                vel(i, j, k, 2) = 0.0;

                const amrex::Real f =
                    std::pow(x, 4) - 2 * std::pow(x, 3) + std::pow(x, 2);
                const amrex::Real f1 =
                    4 * std::pow(x, 3) - 6 * std::pow(x, 2) + 2 * x;
                const amrex::Real f3 = 24 * x - 12;
                const amrex::Real g = std::pow(y, 4) - std::pow(y, 2);
                const amrex::Real g1 = 4 * std::pow(y, 3) - 2 * y;
                const amrex::Real g2 = 12 * std::pow(y, 2) - 2;

                const amrex::Real F = std::pow(x, 5) / 5 - std::pow(x, 4) / 2 +
                                      std::pow(x, 3) / 3;
                const amrex::Real F1 = -4 * std::pow(x, 6) +
                                       12 * std::pow(x, 5) -
                                       14 * std::pow(x, 4) +
                                       8 * std::pow(x, 3) - 2 * std::pow(x, 2);
                const amrex::Real F2 = f * f / 2;
                const amrex::Real G1 =
                    -24 * std::pow(y, 5) + 8 * std::pow(y, 3) - 4 * y;

                src(i, j, k, 0) = 0.0;
                src(i, j, k, 1) =
                    (8 / Re * (24 * F + 2 * f1 * g2 + f3 * g) +
                     64 * (F2 * G1 - g * g1 * F1));
                src(i, j, k, 2) = 0.0;
            });
    }
}

template <typename T>
amrex::Real BurggrafFlow::compute_error(const Field& field)
{
    amrex::Real error = 0.0;
    T f_exact;
    const auto comp = f_exact.m_comp;

    const int nlevels = m_repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {

        amrex::iMultiFab level_mask;
        if (lev < nlevels - 1) {
            level_mask = makeFineMask(
                m_mesh.boxArray(lev), m_mesh.DistributionMap(lev),
                m_mesh.boxArray(lev + 1), m_mesh.refRatio(lev), 1, 0);
        } else {
            level_mask.define(
                m_mesh.boxArray(lev), m_mesh.DistributionMap(lev), 1, 0,
                amrex::MFInfo());
            level_mask.setVal(1);
        }

        const auto& dx = m_mesh.Geom(lev).CellSizeArray();
        const auto& prob_lo = m_mesh.Geom(lev).ProbLoArray();

        const auto& fld = field(lev);
        auto const& fld_arr = fld.const_arrays();
        auto const& mask_arr = level_mask.const_arrays();

        error += amrex::ParReduce(
            amrex::TypeList<amrex::ReduceOpSum>{},
            amrex::TypeList<amrex::Real>{}, fld, amrex::IntVect(0),
            [=] AMREX_GPU_HOST_DEVICE(int box_no, int i, int j, int k)
                -> amrex::GpuTuple<amrex::Real> {
                auto const& fld_bx = fld_arr[box_no];
                auto const& mask_bx = mask_arr[box_no];

                amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
                amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
                const amrex::Real u = fld_bx(i, j, k, comp);
                const amrex::Real u_exact = f_exact(x, y);
                const amrex::Real cell_vol = dx[0] * dx[1] * dx[2];

                return cell_vol * mask_bx(i, j, k) * (u - u_exact) *
                       (u - u_exact);
            });
    }

    amrex::ParallelDescriptor::ReduceRealSum(error);

    const amrex::Real total_vol = m_mesh.Geom(0).ProbDomain().volume();
    return std::sqrt(error / total_vol);
}

void BurggrafFlow::output_error()
{
    // TODO: gradp analytical solution has not been adjusted for mesh
    // mapping
    const amrex::Real u_err = compute_error<UExact>(m_velocity);
    const amrex::Real v_err = compute_error<VExact>(m_velocity);

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f;
        f.open(m_output_fname.c_str(), std::ios_base::app);
        f << std::setprecision(12) << m_time.new_time() << std::setw(m_w)
          << std::setw(m_w) << u_err << std::setw(m_w) << v_err << std::endl;
        f.close();
    }
}

void BurggrafFlow::post_init_actions() { output_error(); }

void BurggrafFlow::post_advance_work() { output_error(); }

} // namespace amr_wind::burggraf
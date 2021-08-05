#include "amr-wind/utilities/sampling/Enstrophy.H"
#include "amr-wind/utilities/io_utils.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include <AMReX_MultiFabUtil.H>
#include "AMReX_ParmParse.H"
#include "amr-wind/utilities/IOManager.H"
#include "amr-wind/fvm/vorticity_mag.H"

namespace amr_wind {
namespace enstrophy {

Enstrophy::Enstrophy(CFDSim& sim, const std::string& label)
    : m_sim(sim)
    , m_label(label)
    , m_velocity(sim.repo().get_field("velocity"))
    , m_density(sim.repo().get_field("density"))
{}

Enstrophy::~Enstrophy() = default;

void Enstrophy::initialize()
{
    BL_PROFILE("amr-wind::Enstrophy::initialize");
    amrex::ParmParse pp(m_label);
    pp.query("output_frequency", m_out_freq);

    prepare_ascii_file();
}

amrex::Real Enstrophy::calculate_enstrophy()
{
    BL_PROFILE("amr-wind::Enstrophy::calculate_enstrophy");

    // integrated total Enstrophy
    amrex::Real total_enstrophy = 0.0;

    const int finest_level = m_velocity.repo().num_active_levels() - 1;
    const auto& geom = m_velocity.repo().mesh().Geom();

    auto vorticity = amr_wind::fvm::vorticity_mag(m_velocity);

    for (int lev = 0; lev <= finest_level; lev++) {

        amrex::iMultiFab level_mask;
        if (lev < finest_level) {
            level_mask = makeFineMask(
                m_sim.mesh().boxArray(lev), m_sim.mesh().DistributionMap(lev),
                m_sim.mesh().boxArray(lev + 1), amrex::IntVect(2), 1, 0);
        } else {
            level_mask.define(
                m_sim.mesh().boxArray(lev), m_sim.mesh().DistributionMap(lev),
                1, 0, amrex::MFInfo());
            level_mask.setVal(1);
        }

        const amrex::Real cell_vol = geom[lev].CellSize()[0] *
                                     geom[lev].CellSize()[1] *
                                     geom[lev].CellSize()[2];

        total_enstrophy += amrex::ReduceSum(
            m_density(lev), (*vorticity)(lev), level_mask, 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& den_arr,
                amrex::Array4<amrex::Real const> const& vort_arr,
                amrex::Array4<int const> const& mask_arr) -> amrex::Real {
                amrex::Real enstrophy_fab = 0.0;

                amrex::Loop(
                    bx, [=, &enstrophy_fab](int i, int j, int k) noexcept {
                        enstrophy_fab +=
                            cell_vol * mask_arr(i, j, k) * den_arr(i, j, k) *
                            (vort_arr(i, j, k) * vort_arr(i, j, k));
                    });
                return enstrophy_fab;
            });
    }

    // total volume of grid on level 0
    const amrex::Real total_vol = geom[0].ProbDomain().volume();

    total_enstrophy *= 0.5 / total_vol;

    amrex::ParallelDescriptor::ReduceRealSum(total_enstrophy);

    return total_enstrophy;
}

void Enstrophy::post_advance_work()
{
    BL_PROFILE("amr-wind::Enstrophy::post_advance_work");
    const auto& time = m_sim.time();
    const int tidx = time.time_index();
    // Skip processing if it is not an output timestep
    if (!(tidx % m_out_freq == 0)) return;

    m_total_enstrophy = calculate_enstrophy();

    write_ascii();
}

void Enstrophy::prepare_ascii_file()
{
    BL_PROFILE("amr-wind::Enstrophy::prepare_ascii_file");

    const std::string post_dir = "post_processing";
    const std::string sname =
        amrex::Concatenate(m_label, m_sim.time().time_index());

    if (!amrex::UtilCreateDirectory(post_dir, 0755)) {
        amrex::CreateDirectoryFailed(post_dir);
    }
    m_out_fname = post_dir + "/" + sname + ".txt";

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f(m_out_fname.c_str());
        f << "time_step time enstrophy" << std::endl;
        f.close();
    }
}

void Enstrophy::write_ascii()
{
    BL_PROFILE("amr-wind::Enstrophy::write_ascii");

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f(m_out_fname.c_str(), std::ios_base::app);
        f << m_sim.time().time_index() << std::scientific
          << std::setprecision(m_precision) << std::setw(m_width)
          << m_sim.time().new_time();
        f << std::setw(m_width) << m_total_enstrophy;
        f << std::endl;
        f.close();
    }
}

} // namespace enstrophy
} // namespace amr_wind

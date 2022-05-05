#include "amr-wind/utilities/sampling/WaveEnergy.H"
#include "amr-wind/utilities/io_utils.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include <AMReX_MultiFabUtil.H>
#include <utility>
#include "AMReX_ParmParse.H"
#include "amr-wind/utilities/IOManager.H"

namespace amr_wind {
namespace wave_energy {

WaveEnergy::WaveEnergy(CFDSim& sim, std::string label)
    : m_sim(sim)
    , m_label(std::move(label))
    , m_velocity(sim.repo().get_field("velocity"))
    , m_vof(sim.repo().get_field("vof"))
{}

WaveEnergy::~WaveEnergy() = default;

void WaveEnergy::initialize()
{
    BL_PROFILE("amr-wind::WaveEnergy::initialize");
    {
        amrex::ParmParse pp(m_label);
        pp.query("output_frequency", m_out_freq);
        pp.query("potential_energy_offset", m_pe_off);
    }
    // Get gravity and density
    {
        amrex::ParmParse pp("incflo");
        pp.queryarr("gravity", m_gravity);
    }
    auto& mphase = m_sim.physics_manager().get<amr_wind::MultiPhase>();
    m_rho1 = mphase.rho1();

    prepare_ascii_file();
}

amrex::Real WaveEnergy::calculate_kinetic_energy()
{
    BL_PROFILE("amr-wind::WaveEnergy::calculate_kinetic_energy");

    // integrated total wave Energy
    amrex::Real wave_ke = 0.0;

    // Liquid density
    const amrex::Real rho_liq = m_rho1;

    const int finest_level = m_velocity.repo().num_active_levels() - 1;
    const auto& geom = m_velocity.repo().mesh().Geom();

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

        wave_ke += amrex::ReduceSum(
            m_vof(lev), m_velocity(lev), level_mask, 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& vof_arr,
                amrex::Array4<amrex::Real const> const& vel_arr,
                amrex::Array4<int const> const& mask_arr) -> amrex::Real {
                amrex::Real Wave_Energy_Fab = 0.0;

                amrex::Loop(
                    bx, [=, &Wave_Energy_Fab](int i, int j, int k) noexcept {
                        Wave_Energy_Fab +=
                            cell_vol * mask_arr(i, j, k) * 0.5 * rho_liq *
                            vof_arr(i, j, k) *
                            (vel_arr(i, j, k, 0) * vel_arr(i, j, k, 0) +
                             vel_arr(i, j, k, 1) * vel_arr(i, j, k, 1) +
                             vel_arr(i, j, k, 2) * vel_arr(i, j, k, 2));
                    });
                return Wave_Energy_Fab;
            });
    }

    amrex::ParallelDescriptor::ReduceRealSum(wave_ke);

    return wave_ke;
}

amrex::Real WaveEnergy::calculate_potential_energy()
{
    BL_PROFILE("amr-wind::WaveEnergy::calculate_potential_energy");

    // integrated total wave Energy
    amrex::Real wave_pe = 0.0;
    // Liquid density
    const amrex::Real rho_liq = m_rho1;
    // gravity constant
    const amrex::Real g = -m_gravity[2];

    const int finest_level = m_velocity.repo().num_active_levels() - 1;
    const auto& geom = m_velocity.repo().mesh().Geom();

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

        const amrex::Real dz = geom[lev].CellSize()[2];
        const amrex::Real probloz = geom[lev].ProbLo()[2];

        wave_pe += amrex::ReduceSum(
            m_vof(lev), level_mask, 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& vof_arr,
                amrex::Array4<int const> const& mask_arr) -> amrex::Real {
                amrex::Real Wave_Energy_Fab = 0.0;

                amrex::Loop(
                    bx, [=, &Wave_Energy_Fab](int i, int j, int k) noexcept {
                        // Crude model of liquid height in multiphase cells
                        amrex::Real kk =
                            (vof_arr(i, j, k + 1) > vof_arr(i, j, k)) ? k + 1
                                                                      : k;
                        amrex::Real dir =
                            (vof_arr(i, j, k + 1) > vof_arr(i, j, k)) ? -1 : 1;
                        const amrex::Real zl =
                            probloz + (kk + dir * 0.5 * vof_arr(i, j, k)) * dz;
                        Wave_Energy_Fab += cell_vol * mask_arr(i, j, k) *
                                           rho_liq * vof_arr(i, j, k) * g * zl;
                    });
                return Wave_Energy_Fab;
            });
    }

    amrex::ParallelDescriptor::ReduceRealSum(wave_pe);

    return wave_pe;
}

void WaveEnergy::post_advance_work()
{
    BL_PROFILE("amr-wind::WaveEnergy::post_advance_work");
    const auto& time = m_sim.time();
    const int tidx = time.time_index();
    // Skip processing if it is not an output timestep
    if (!(tidx % m_out_freq == 0)) {
        return;
    }

    m_wave_kinetic_energy = calculate_kinetic_energy();
    m_wave_potential_energy = calculate_potential_energy() + m_pe_off;

    write_ascii();
}

void WaveEnergy::prepare_ascii_file()
{
    BL_PROFILE("amr-wind::WaveEnergy::prepare_ascii_file");

    const std::string post_dir = "post_processing";
    const std::string sname =
        amrex::Concatenate(m_label, m_sim.time().time_index());

    if (!amrex::UtilCreateDirectory(post_dir, 0755)) {
        amrex::CreateDirectoryFailed(post_dir);
    }
    m_out_fname = post_dir + "/" + sname + ".txt";

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f(m_out_fname.c_str());
        f << "time_step time kinetic_energy potential_energy" << std::endl;
        f.close();
    }
}

void WaveEnergy::write_ascii()
{
    BL_PROFILE("amr-wind::WaveEnergy::write_ascii");

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f(m_out_fname.c_str(), std::ios_base::app);
        f << m_sim.time().time_index() << std::scientific
          << std::setprecision(m_precision) << std::setw(m_width)
          << m_sim.time().new_time();
        f << std::setw(m_width) << m_wave_kinetic_energy;
        f << std::setw(m_width) << m_wave_potential_energy;
        f << std::endl;
        f.close();
    }
}

} // namespace wave_energy
} // namespace amr_wind

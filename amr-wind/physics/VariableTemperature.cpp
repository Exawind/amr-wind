#include "amr-wind/physics/VariableTemperature.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParmParse.H"
#include "AMReX_ParReduce.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/utilities/IOManager.H"

namespace amr_wind::VariableTemperature {

namespace {} // namespace

VariableTemperature::VariableTemperature(CFDSim& sim)
    : m_time(sim.time())
    , m_sim(sim)
    , m_repo(sim.repo())
    , m_mesh(sim.mesh())
    , m_velocity(sim.repo().get_field("velocity"))
    , m_surf_temp(sim.repo().declare_field("surf_temp", 3, 1, 1))

{
    amrex::ParmParse pp("VariableTemperature");
    pp.query("data_spacing", m_separation);
    m_current_time = m_time.new_time();
    // Read the first file
    std::string groundfile("ground_mesh.txt");
    std::ifstream file1(groundfile, std::ios::in);
    if (!file1.good()) {
        amrex::Abort("Cannot find file ground_mesh.txt");
    }
    amrex::Real value1, value2, value3;
    while (file1 >> value1 >> value2) {
        m_xloc.push_back(value1);
        m_yloc.push_back(value2);
    }
    file1.close();
    device_xloc.resize(m_xloc.size());
    device_yloc.resize(m_yloc.size());
    std::string temperaturefile("ground_temperature.txt");
    std::ifstream file2(temperaturefile, std::ios::in);
    if (!file2.good()) {
        amrex::Abort("Cannot find file ground_temperature.txt");
    }
    while (file2 >> value3) {
        m_Tloc.push_back(value3);
    }
    device_Tloc.resize(m_Tloc.size());
    m_sim.io_manager().register_io_var("surf_temp");
}

void VariableTemperature::post_init_actions()
{
    BL_PROFILE("amr-wind::" + this->identifier() + "::post_init_actions");
    const auto& geom_vec = m_sim.repo().mesh().Geom();
    const int nlevels = m_sim.repo().num_active_levels();
    for (int level = 0; level < nlevels; ++level) {
        const auto& geom = geom_vec[level];
        const auto& dx = geom.CellSizeArray();
        const auto& prob_lo = geom.ProbLoArray();
        auto& velocity = m_velocity(level);
        auto& surf_temp = m_surf_temp(level);
        // Copy Temperature to gpu
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_xloc.begin(), m_xloc.end(),
            device_xloc.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_yloc.begin(), m_yloc.end(),
            device_yloc.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_Tloc.begin(), m_Tloc.end(),
            device_Tloc.begin());
        const auto* xloc_ptr = device_xloc.data();
        const auto* yloc_ptr = device_yloc.data();
        const auto* Tloc_ptr = device_Tloc.data();
        for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
            const auto& vbx = mfi.validbox();
            auto level_surf_temp = surf_temp.array(mfi);
            const unsigned mesh_size = m_xloc.size();
            const int index =
                std::max(0, int(m_current_time / m_separation) - 1);
            amrex::ParallelFor(
                vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    // compute the source term
                    const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
                    const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
                    amrex::Real residual = 10000;
                    amrex::Real locT1 = 300;
                    amrex::Real locT2 = 300;
                    for (unsigned ii = 0; ii < mesh_size; ++ii) {
                        const amrex::Real radius = std::sqrt(
                            std::pow(x - xloc_ptr[ii], 2) +
                            std::pow(y - yloc_ptr[ii], 2));
                        if (radius < residual) {
                            residual = radius;
                            locT1 = Tloc_ptr[index * mesh_size + ii];
                            locT2 = Tloc_ptr[(index + 1) * mesh_size + ii];
                        }
                    }
                    level_surf_temp(i, j, k, 0) = locT1; // Lower Value
                    level_surf_temp(i, j, k, 1) = locT1; // Actual Used Value
                    level_surf_temp(i, j, k, 2) = locT2; // Upper Value
                });
        }
    }
}

void VariableTemperature::pre_init_actions()
{
    BL_PROFILE("amr-wind::" + this->identifier() + "::pre_init_actions");
}

void VariableTemperature::pre_advance_work()
{
    BL_PROFILE("amr-wind::" + this->identifier() + "::pre_advanced_work");
    const amrex::Real fraction = std::max(
        0.0,
        std::min((m_current_time - m_time.new_time()) / m_separation, 1.0));
    if (m_current_time >= m_time.new_time()) {
        const int nlevels = m_sim.repo().num_active_levels();
        for (int level = 0; level < nlevels; ++level) {
            auto& velocity = m_velocity(level);
            auto& surf_temp = m_surf_temp(level);
            for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
                const auto& vbx = mfi.validbox();
                auto level_surf_temp = surf_temp.array(mfi);
                amrex::ParallelFor(
                    vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        level_surf_temp(i, j, k, 1) =
                            fraction * level_surf_temp(i, j, k, 0) +
                            (1 - fraction) * level_surf_temp(i, j, k, 2);
                    });
            }
        }
        return;
    }
    m_current_time += m_separation;
    const auto& geom_vec = m_sim.repo().mesh().Geom();
    const int nlevels = m_sim.repo().num_active_levels();
    for (int level = 0; level < nlevels; ++level) {
        const auto& geom = geom_vec[level];
        const auto& dx = geom.CellSizeArray();
        const auto& prob_lo = geom.ProbLoArray();
        auto& velocity = m_velocity(level);
        auto& surf_temp = m_surf_temp(level);
        const auto* xloc_ptr = device_xloc.data();
        const auto* yloc_ptr = device_yloc.data();
        const auto* Tloc_ptr = device_Tloc.data();
        for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
            const auto& vbx = mfi.validbox();
            auto level_surf_temp = surf_temp.array(mfi);
            const unsigned mesh_size = m_xloc.size();
            const int index =
                std::max(0, int(m_current_time / m_separation) - 1);
            amrex::ParallelFor(
                vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
                    const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
                    amrex::Real residual = 10000;
                    amrex::Real locT1 = 300;
                    amrex::Real locT2 = 300;
                    for (unsigned ii = 0; ii < mesh_size; ++ii) {
                        const amrex::Real radius = std::sqrt(
                            std::pow(x - xloc_ptr[ii], 2) +
                            std::pow(y - yloc_ptr[ii], 2));
                        if (radius < residual) {
                            residual = radius;
                            locT1 = Tloc_ptr[index * mesh_size + ii];
                            locT2 = Tloc_ptr[(index + 1) * mesh_size + ii];
                        }
                    }
                    level_surf_temp(i, j, k, 0) = locT1;
                    level_surf_temp(i, j, k, 2) = locT2;
                    level_surf_temp(i, j, k, 1) =
                        fraction * level_surf_temp(i, j, k, 0) +
                        (1 - fraction) * level_surf_temp(i, j, k, 2);
                });
        }
    }
}

} // namespace amr_wind::VariableTemperature

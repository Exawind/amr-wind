#include "amr-wind/physics/MultiTemperature.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParmParse.H"
#include "AMReX_ParReduce.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/utilities/IOManager.H"

namespace amr_wind::multitemperature {

namespace {} // namespace

MultiTemperature::MultiTemperature(CFDSim& sim)
    : m_time(sim.time())
    , m_sim(sim)
    , m_repo(sim.repo())
    , m_mesh(sim.mesh())
    , m_velocity(sim.repo().get_field("velocity"))
    , m_surf_temp(sim.repo().declare_field("surf_temp", 1, 1, 1))
    , m_surf_temp_lower(sim.repo().declare_field("surf_temp_lower", 1, 1, 1))
    , m_surf_temp_upper(sim.repo().declare_field("surf_temp_upper", 1, 1, 1))
{
    amrex::ParmParse pp("MultiTemperature");
    pp.query("data_spacing", m_separation);
    m_current_time = m_time.new_time();
    // Read the first file
    std::string groundfile("ground.mesh");
    std::ifstream file1(groundfile, std::ios::in);
    if (!file1.good()) {
        amrex::Abort("Cannot find ground.meshfile");
    }
    amrex::Real value1, value2, value3;
    while (file1 >> value1 >> value2) {
        m_xloc.push_back(value1);
        m_yloc.push_back(value2);
    }
    file1.close();
    device_xloc.resize(m_xloc.size());
    device_yloc.resize(m_yloc.size());
    std::string temperaturefile("ground.temperature");
    std::ifstream file2(temperaturefile, std::ios::in);
    if (!file1.good()) {
        amrex::Abort("Cannot find ground.temperature");
    }
    while (file2 >> value3) {
        m_Tloc.push_back(value3);
    }
    device_Tloc.resize(m_Tloc.size());
    m_sim.io_manager().register_io_var("surf_temp");
}

void MultiTemperature::post_init_actions()
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
        auto& surf_temp_lower = m_surf_temp_lower(level);
        auto& surf_temp_upper = m_surf_temp_upper(level);
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
            auto level_surf_temp_lower = surf_temp_lower.array(mfi);
            auto level_surf_temp_upper = surf_temp_upper.array(mfi);
            const unsigned mesh_size = m_xloc.size();
            const int index =
                std::max(0, int(m_current_time / m_separation) - 1);
            amrex::Print() << "Checking issues:" << m_Tloc.size() << "  "
                           << mesh_size << "   " << index
                           << index * (mesh_size - 1) << "  "
                           << (index + 1) * (mesh_size - 1) << std::endl;
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
                            locT1 = Tloc_ptr[index * (mesh_size - 1) + ii];
                            locT2 =
                                Tloc_ptr[(index + 1) * (mesh_size - 1) + ii];
                        }
                    }
                    level_surf_temp(i, j, k, 0) = locT1;
                    level_surf_temp_lower(i, j, k, 0) = locT1;
                    level_surf_temp_upper(i, j, k, 0) = locT2;
                });
        }
    }
}

void MultiTemperature::pre_init_actions()
{
    BL_PROFILE("amr-wind::" + this->identifier() + "::pre_init_actions");
}

void MultiTemperature::pre_advance_work()
{
    BL_PROFILE("amr-wind::" + this->identifier() + "::pre_advanced_work");
    if (m_current_time >= m_time.new_time()) {
        const amrex::Real fraction = std::max(
            0.0,
            std::min((m_current_time - m_time.new_time()) / m_separation, 1.0));
        amrex::Print() << "Interpolation:" << m_time.new_time() << "  "
                       << m_current_time << "   " << fraction << std::endl;
        const int nlevels = m_sim.repo().num_active_levels();
        for (int level = 0; level < nlevels; ++level) {
            auto& velocity = m_velocity(level);
            auto& surf_temp = m_surf_temp(level);
            auto& surf_temp_lower = m_surf_temp_lower(level);
            auto& surf_temp_upper = m_surf_temp_upper(level);
            for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
                const auto& vbx = mfi.validbox();
                auto level_surf_temp = surf_temp.array(mfi);
                auto level_surf_temp_lower = surf_temp_lower.array(mfi);
                auto level_surf_temp_upper = surf_temp_upper.array(mfi);
                amrex::ParallelFor(
                    vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        level_surf_temp(i, j, k) =
                            fraction * level_surf_temp_lower(i, j, k) +
                            (1 - fraction) * level_surf_temp_upper(i, j, k);
                    });
            }
        }
        return;
    }
    amrex::Print() << "Changing:" << std::endl;
    m_current_time += m_separation;
    const auto& geom_vec = m_sim.repo().mesh().Geom();
    const int nlevels = m_sim.repo().num_active_levels();
    for (int level = 0; level < nlevels; ++level) {
        const auto& geom = geom_vec[level];
        const auto& dx = geom.CellSizeArray();
        const auto& prob_lo = geom.ProbLoArray();
        auto& velocity = m_velocity(level);
        auto& surf_temp = m_surf_temp(level);
        auto& surf_temp_lower = m_surf_temp_lower(level);
        auto& surf_temp_upper = m_surf_temp_upper(level);
        const auto* xloc_ptr = device_xloc.data();
        const auto* yloc_ptr = device_yloc.data();
        const auto* Tloc_ptr = device_Tloc.data();
        for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
            const auto& vbx = mfi.validbox();
            auto level_surf_temp = surf_temp.array(mfi);
            auto level_surf_temp_lower = surf_temp_lower.array(mfi);
            auto level_surf_temp_upper = surf_temp_upper.array(mfi);
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
                            locT1 = Tloc_ptr[index * (mesh_size - 1) + ii];
                            locT2 =
                                Tloc_ptr[(index + 1) * (mesh_size - 1) + ii];
                        }
                    }
                    level_surf_temp(i, j, k, 0) = locT1;
                    level_surf_temp_lower(i, j, k, 0) = locT1;
                    level_surf_temp_upper(i, j, k, 0) = locT2;
                });
        }
    }
}

} // namespace amr_wind::multitemperature

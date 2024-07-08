#include "amr-wind/physics/TerrainDrag.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParmParse.H"
#include "AMReX_ParReduce.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/utilities/IOManager.H"

namespace amr_wind::terraindrag {

namespace {} // namespace

TerrainDrag::TerrainDrag(CFDSim& sim)
    : m_sim(sim)
    , m_repo(sim.repo())
    , m_mesh(sim.mesh())
    , m_velocity(sim.repo().get_field("velocity"))
    , m_terrain_blank(sim.repo().declare_field("terrainBlank", 1, 1, 1))
    , m_terrain_drag(sim.repo().declare_field("terrainDrag", 1, 1, 1))
    , m_terrainz0(sim.repo().declare_field("terrainz0", 1, 1, 1))
{
    std::string terrainfile("terrain.amrwind");
    std::ifstream file(terrainfile, std::ios::in);
    if (!file.good()) {
        amrex::Abort("Cannot find terrain.amrwind  file");
    }
    amrex::Real value1, value2, value3;
    while (file >> value1 >> value2 >> value3) {
        m_xterrain.push_back(value1);
        m_yterrain.push_back(value2);
        m_zterrain.push_back(value3);
    }
    file.close();
    // No checks for the file as it is optional currently
    std::string roughnessfile("terrain.roughness");
    std::ifstream file1(roughnessfile, std::ios::in);
    while (file1 >> value1 >> value2 >> value3) {
        m_xrough.push_back(value1);
        m_yrough.push_back(value2);
        m_z0rough.push_back(value3);
    }
    file1.close();
    m_sim.io_manager().register_io_var("terrainDrag");
    m_sim.io_manager().register_io_var("terrainBlank");
    m_sim.io_manager().register_io_var("terrainz0");
}

void TerrainDrag::post_init_actions()
{
    BL_PROFILE("amr-wind::" + this->identifier() + "::post_init_actions");
    const auto& geom_vec = m_sim.repo().mesh().Geom();
    const int nlevels = m_sim.repo().num_active_levels();
    for (int level = 0; level < nlevels; ++level) {
        const auto& geom = geom_vec[level];
        const auto& dx = geom.CellSizeArray();
        const auto& prob_lo = geom.ProbLoArray();
        auto& velocity = m_velocity(level);
        auto& blanking = m_terrain_blank(level);
        auto& terrainz0 = m_terrainz0(level);
        auto& drag = m_terrain_drag(level);
        // copy terrain data to gpu
        amrex::Gpu::DeviceVector<amrex::Real> device_xterrain(m_xterrain.size());
        amrex::Gpu::DeviceVector<amrex::Real> device_yterrain(m_xterrain.size());
        amrex::Gpu::DeviceVector<amrex::Real> device_zterrain(m_xterrain.size());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_xterrain.begin(), m_xterrain.end(),
            device_xterrain.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_yterrain.begin(), m_yterrain.end(),
            device_yterrain.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_zterrain.begin(), m_zterrain.end(),
            device_zterrain.begin());
        const auto* xterrain_ptr = device_xterrain.data();
        const auto* yterrain_ptr = device_yterrain.data();
        const auto* zterrain_ptr = device_zterrain.data();
        // Copy Roughness to gpu
        amrex::Gpu::DeviceVector<amrex::Real> device_xrough(m_xrough.size());
        amrex::Gpu::DeviceVector<amrex::Real> device_yrough(m_xrough.size());
        amrex::Gpu::DeviceVector<amrex::Real> device_z0rough(m_xrough.size());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_xrough.begin(), m_xrough.end(),
            device_xrough.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_yrough.begin(), m_yrough.end(),
            device_yrough.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_z0rough.begin(), m_z0rough.end(),
            device_z0rough.begin());
        const auto* xrough_ptr = device_xrough.data();
        const auto* yrough_ptr = device_yrough.data();
        const auto* z0rough_ptr = device_z0rough.data();
        for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
            const auto& vbx = mfi.validbox();
            auto levelBlanking = blanking.array(mfi);
            auto levelDrag = drag.array(mfi);
            auto levelz0 = terrainz0.array(mfi);
            const unsigned terrainSize = m_xterrain.size();
            const unsigned roughnessSize = m_xrough.size();
            amrex::ParallelFor(
                vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    // compute the source term
                    const amrex::Real x1 = prob_lo[0] + (i + 0.5) * dx[0];
                    const amrex::Real x2 = prob_lo[1] + (j + 0.5) * dx[1];
                    const amrex::Real x3 = prob_lo[2] + (k + 0.5) * dx[2];
                    // Terrain Height
                    amrex::Real residual = 10000;
                    amrex::Real terrainHt = 0.0;
                    for (unsigned ii = 0; ii < terrainSize; ++ii) {
                        const amrex::Real radius = std::sqrt(
                            std::pow(x1 - xterrain_ptr[ii], 2) +
                            std::pow(x2 - yterrain_ptr[ii], 2));
                        if (radius < residual) {
                            residual = radius;
                            terrainHt = zterrain_ptr[ii];
                        }
                        if (radius < dx[0]) {
                            break;
                        }
                    }
                    levelBlanking(i, j, k, 0) = static_cast<float>(x3 <= terrainHt);
                    residual = 10000;
                    amrex::Real roughz0 = 0.1;
                    for (unsigned ii = 0; ii < roughnessSize; ++ii) {
                        const amrex::Real radius = std::sqrt(
                            std::pow(x1 - xrough_ptr[ii], 2) +
                            std::pow(x2 - yrough_ptr[ii], 2));
                        if (radius < residual) {
                            residual = radius;
                            roughz0 = z0rough_ptr[ii];
                        }
                        if (radius < dx[0]) {
                            break;
                        }
                    }
                    levelz0(i, j, k, 0) = roughz0;
                });
            amrex::ParallelFor(
                vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    // Terrain Height
                    amrex::Real residual = 10000;
                    amrex::Real terrainHt = 0.0;
                    const amrex::Real x1 = prob_lo[0] + (i + 0.5) * dx[0];
                    const amrex::Real x2 = prob_lo[1] + (j + 0.5) * dx[1];
                    const amrex::Real x3 = prob_lo[2] + (k + 0.5) * dx[2];
                    for (unsigned ii = 0; ii < terrainSize; ++ii) {
                        const amrex::Real radius = std::sqrt(
                            std::pow(x1 - xterrain_ptr[ii], 2) +
                            std::pow(x2 - yterrain_ptr[ii], 2));
                        if (radius < residual) {
                            residual = radius;
                            terrainHt = zterrain_ptr[ii];
                        }
                    }
                    levelDrag(i, j, k, 0) = 0;
                    if (x3 > terrainHt && k > 0 &&
                        levelBlanking(i, j, k - 1, 0) == 1) {
                        levelDrag(i, j, k, 0) = 1;
                    }
                });
        }
    }
}

void TerrainDrag::pre_init_actions()
{
    BL_PROFILE("amr-wind::" + this->identifier() + "::pre_init_actions");
}

int TerrainDrag::returnBlankValue(int i, int j, int k)
{
    int lev = 0;
    amrex::MFIter mfi(m_terrain_blank(lev));
    const auto& levelBlanking = m_terrain_blank(lev).const_array(mfi);
    return int(levelBlanking(i, j, k));
}

} // namespace amr_wind::terraindrag

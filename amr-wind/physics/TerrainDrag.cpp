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
    , m_terrainBlank(sim.repo().declare_field("terrainBlank", 1, 1, 1))
    , m_terrainDrag(sim.repo().declare_field("terrainDrag", 1, 1, 1))
{
    std::ifstream file("terrain.amrwind");
    amrex::Real value1, value2, value3;
    while (file >> value1 >> value2 >> value3) {
        m_xterrain.push_back(value1);
        m_yterrain.push_back(value2);
        m_zterrain.push_back(value3);
    }
    file.close();
    m_sim.io_manager().register_io_var("terrainDrag");
    m_sim.io_manager().register_io_var("terrainBlank");
}

void TerrainDrag::initialize_fields(int level, const amrex::Geometry& geom)
{
    using namespace utils;

    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();
    auto& velocity = m_velocity(level);
    auto& blanking = m_terrainBlank(level);
    auto& drag = m_terrainDrag(level);
    // amrex::Print()<<" Terrain Initial:"<<level<<std::endl;
    //  copy terrain data to gpu
    amrex::Gpu::DeviceVector<amrex::Real> gpu_xterrain(m_xterrain.size());
    amrex::Gpu::DeviceVector<amrex::Real> gpu_yterrain(m_xterrain.size());
    amrex::Gpu::DeviceVector<amrex::Real> gpu_zterrain(m_xterrain.size());
    auto* xterrain_ptr = gpu_xterrain.data();
    auto* yterrain_ptr = gpu_yterrain.data();
    auto* zterrain_ptr = gpu_zterrain.data();
    amrex::ParallelFor(m_xterrain.size(), [=] AMREX_GPU_DEVICE(int n) {
        xterrain_ptr[n] = m_xterrain[n];
        yterrain_ptr[n] = m_yterrain[n];
        zterrain_ptr[n] = m_zterrain[n];
    });
    for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        auto levelBlanking = blanking.array(mfi);
        auto levelDrag = drag.array(mfi);
        const amrex::Real terrainSize = m_xterrain.size();
        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // compute the source term
                const amrex::Real x1 = prob_lo[0] + (i + 0.5) * dx[0];
                const amrex::Real x2 = prob_lo[1] + (j + 0.5) * dx[1];
                const amrex::Real x3 = prob_lo[2] + (k + 0.5) * dx[2];
                // Terrain Height
                amrex::Real residual = 10000;
                amrex::Real terrainHt = 0.0;
                for (int ii = 0; ii < terrainSize; ++ii) {
                    const amrex::Real radius = std::sqrt(
                        std::pow(x1 - xterrain_ptr[ii], 2) +
                        std::pow(x2 - yterrain_ptr[ii], 2));
                    if (radius < residual) {
                        residual = radius;
                        terrainHt = zterrain_ptr[ii];
                    }
                }
                const amrex::Real turnOn = (x3 <= terrainHt) ? 1.0 : 0.0;
                levelBlanking(i, j, k, 0) = turnOn;
            });
        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // Terrain Height
                amrex::Real residual = 10000;
                amrex::Real terrainHt = 0.0;
                const amrex::Real x1 = prob_lo[0] + (i + 0.5) * dx[0];
                const amrex::Real x2 = prob_lo[1] + (j + 0.5) * dx[1];
                const amrex::Real x3 = prob_lo[2] + (k + 0.5) * dx[2];
                for (int ii = 0; ii < terrainSize; ++ii) {
                    const amrex::Real radius = std::sqrt(
                        std::pow(x1 - xterrain_ptr[ii], 2) +
                        std::pow(x2 - yterrain_ptr[ii], 2));
                    if (radius < residual) {
                        residual = radius;
                        terrainHt = zterrain_ptr[ii];
                    }
                }
                levelDrag(i, j, k, 0) = 0.0;
                if (x3 > terrainHt && k > 0 &&
                    levelBlanking(i, j, k - 1, 0) == 1) {
                    //	    amrex::Print()<<"Adding Drag:"<<x1<<"  "<<x2<<"
                    //"<<x3<<"   "<<levelBlanking(i,j,k-1,0)<<std::endl;
                    levelDrag(i, j, k, 0) = 1.0;
                }
            });
    }
}

} // namespace amr_wind::terraindrag

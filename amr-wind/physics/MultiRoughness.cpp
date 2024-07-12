#include "amr-wind/physics/MultiRoughness.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParmParse.H"
#include "AMReX_ParReduce.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/utilities/IOManager.H"

namespace amr_wind::multiroughness {

namespace {} // namespace

MultiRoughness::MultiRoughness(CFDSim& sim)
    : m_sim(sim)
    , m_repo(sim.repo())
    , m_mesh(sim.mesh())
    , m_velocity(sim.repo().get_field("velocity"))
    , m_terrainz0(sim.repo().declare_field("lowerz0", 1, 1, 1))
{
    std::string roughnessfile("ground.roughness");
    std::ifstream file1(roughnessfile, std::ios::in);
    if (!file1.good()) {
        amrex::Abort("Cannot find ground.roughness");
    }
    amrex::Real value1, value2, value3;
    while (file1 >> value1 >> value2 >> value3) {
        m_xrough.push_back(value1);
        m_yrough.push_back(value2);
        m_z0rough.push_back(value3);
    }
    file1.close();
    m_sim.io_manager().register_io_var("lowerz0");
}

void MultiRoughness::post_init_actions()
{
    BL_PROFILE("amr-wind::" + this->identifier() + "::post_init_actions");
    const auto& geom_vec = m_sim.repo().mesh().Geom();
    const int nlevels = m_sim.repo().num_active_levels();
    for (int level = 0; level < nlevels; ++level) {
        const auto& geom = geom_vec[level];
        const auto& dx = geom.CellSizeArray();
        const auto& prob_lo = geom.ProbLoArray();
        auto& velocity = m_velocity(level);
        auto& terrainz0 = m_terrainz0(level);
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
            auto levelz0 = terrainz0.array(mfi);
            const unsigned roughnessSize = m_xrough.size();
            amrex::ParallelFor(
                vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    // compute the source term
                    const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
                    const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
                    amrex::Real residual = 10000;
                    amrex::Real roughz0 = 0.1;
                    for (unsigned ii = 0; ii < roughnessSize; ++ii) {
                        const amrex::Real radius = std::sqrt(
                            std::pow(x - xrough_ptr[ii], 2) +
                            std::pow(y - yrough_ptr[ii], 2));
                        if (radius < residual) {
                            residual = radius;
                            roughz0 = z0rough_ptr[ii];
                        }
                    }
                    levelz0(i, j, k, 0) = roughz0;
                });
        }
    }
}

void MultiRoughness::pre_init_actions()
{
    BL_PROFILE("amr-wind::" + this->identifier() + "::pre_init_actions");
}

} // namespace amr_wind::multiroughness

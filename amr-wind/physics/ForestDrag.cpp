#include "amr-wind/physics/ForestDrag.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParmParse.H"
#include "AMReX_ParReduce.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/utilities/IOManager.H"

namespace amr_wind::Forestdrag {

namespace {} // namespace

ForestDrag::ForestDrag(CFDSim& sim)
    : m_sim(sim)
    , m_repo(sim.repo())
    , m_mesh(sim.mesh())
    , m_velocity(sim.repo().get_field("velocity"))
    , m_forest_drag(sim.repo().declare_field("forest_drag", 1, 1, 1))
{
    std::string forestfile("forest.amrwind");
    std::ifstream file(forestfile, std::ios::in);
    if (!file.good()) {
        amrex::Abort("Cannot find forest.amrwind file");
    }
    //! TreeType xc yc height diameter cd lai laimax
    amrex::Real value1, value2, value3, value4, value5, value6, value7, value8;
    while (file >> value1 >> value2 >> value3 >> value4 >> value5 >> value6 >>
           value7 >> value8) {
        m_type_forest.push_back(value1);
        m_x_forest.push_back(value2);
        m_y_forest.push_back(value3);
        m_height_forest.push_back(value4);
        m_diameter_forest.push_back(value5);
        m_cd_forest.push_back(value6);
        m_lai_forest.push_back(value7);
        m_laimax_forest.push_back(value8);
    }
    file.close();
    m_sim.io_manager().register_output_var("forest_drag");
}

void ForestDrag::post_init_actions()
{
    BL_PROFILE("amr-wind::" + this->identifier() + "::post_init_actions");
    const auto& geom_vec = m_sim.repo().mesh().Geom();
    const int nlevels = m_sim.repo().num_active_levels();
    for (int level = 0; level < nlevels; ++level) {
        const auto& geom = geom_vec[level];
        const auto& dx = geom.CellSizeArray();
        const auto& prob_lo = geom.ProbLoArray();
        auto& velocity = m_velocity(level);
        auto& drag = m_forest_drag(level);
        amrex::Gpu::DeviceVector<amrex::Real> device_type_forest(
            m_type_forest.size());
        amrex::Gpu::DeviceVector<amrex::Real> device_x_forest(
            m_x_forest.size());
        amrex::Gpu::DeviceVector<amrex::Real> device_y_forest(
            m_y_forest.size());
        amrex::Gpu::DeviceVector<amrex::Real> device_height_forest(
            m_height_forest.size());
        amrex::Gpu::DeviceVector<amrex::Real> device_diameter_forest(
            m_diameter_forest.size());
        amrex::Gpu::DeviceVector<amrex::Real> device_cd_forest(
            m_cd_forest.size());
        amrex::Gpu::DeviceVector<amrex::Real> device_lai_forest(
            m_lai_forest.size());
        amrex::Gpu::DeviceVector<amrex::Real> device_laimax_forest(
            m_laimax_forest.size());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_type_forest.begin(),
            m_type_forest.end(), device_type_forest.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_x_forest.begin(), m_x_forest.end(),
            device_x_forest.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_y_forest.begin(), m_y_forest.end(),
            device_y_forest.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_height_forest.begin(),
            m_height_forest.end(), device_height_forest.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_diameter_forest.begin(),
            m_diameter_forest.end(), device_diameter_forest.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_cd_forest.begin(), m_cd_forest.end(),
            device_cd_forest.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_lai_forest.begin(), m_lai_forest.end(),
            device_lai_forest.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_laimax_forest.begin(),
            m_laimax_forest.end(), device_laimax_forest.begin());
        const auto* type_forest_ptr = device_type_forest.data();
        const auto* x_forest_ptr = device_x_forest.data();
        const auto* y_forest_ptr = device_y_forest.data();
        const auto* height_forest_ptr = device_height_forest.data();
        const auto* diameter_forest_ptr = device_diameter_forest.data();
        const auto* cd_forest_ptr = device_cd_forest.data();
        const auto* lai_forest_ptr = device_lai_forest.data();
        const auto* laimax_forest_ptr = device_laimax_forest.data();
        for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
            const auto& vbx = mfi.validbox();
            auto levelDrag = drag.array(mfi);
            const auto forestSize = m_x_forest.size();
            amrex::ParallelFor(
                vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    // compute the source term
                    const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
                    const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
                    const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
                    levelDrag(i, j, k) = 0.0;
                    for (unsigned ii = 0; ii < forestSize; ++ii) {
                        const amrex::Real radius = std::sqrt(
                            (x - x_forest_ptr[ii]) * (x - x_forest_ptr[ii]) +
                            (y - y_forest_ptr[ii]) * (y - y_forest_ptr[ii]));
                        amrex::Real af = 0.0;
                        if (z <= height_forest_ptr[ii] &&
                            radius <= (0.5 * diameter_forest_ptr[ii])) {
                            if (type_forest_ptr[ii] == 1) {
                                af = lai_forest_ptr[ii] / height_forest_ptr[ii];
                            } else if (type_forest_ptr[ii] == 2) {
                                amrex::Real ztree = 0;
                                amrex::Real expFun = 0;
                                amrex::Real ratio = 0;
                                const amrex::Real treeZm =
                                    laimax_forest_ptr[ii] *
                                    height_forest_ptr[ii];
                                const amrex::Real dz =
                                    height_forest_ptr[ii] / 100;
                                while (ztree <= height_forest_ptr[ii]) {
                                    ratio = (height_forest_ptr[ii] - treeZm) /
                                            (height_forest_ptr[ii] - ztree);
                                    if (z < treeZm) {
                                        expFun = expFun +
                                                 std::pow(ratio, 6.0) *
                                                     std::exp(6 * (1 - ratio));
                                    } else {
                                        expFun =
                                            expFun +
                                            std::pow(ratio, 0.5) *
                                                std::exp(0.5 * (1 - ratio));
                                    }
                                    ztree = ztree + dz;
                                }
                                const amrex::Real treelaimax =
                                    lai_forest_ptr[ii] / (expFun * dz);
                                ratio = (height_forest_ptr[ii] - treeZm) /
                                        (height_forest_ptr[ii] - z);
                                if (z < treeZm) {
                                    af = treelaimax * std::pow(ratio, 6.0) *
                                         std::exp(6 * (1 - ratio));
                                } else if (z <= height_forest_ptr[ii]) {
                                    af = treelaimax * std::pow(ratio, 0.5) *
                                         std::exp(0.5 * (1 - ratio));
                                }
                            }
                            levelDrag(i, j, k) = cd_forest_ptr[ii] * af;
                        }
                    }
                });
        }
    }
}

void ForestDrag::pre_init_actions()
{
    BL_PROFILE("amr-wind::" + this->identifier() + "::pre_init_actions");
}
} // namespace amr_wind::Forestdrag

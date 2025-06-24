#include "amr-wind/physics/MetMastMarking.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParmParse.H"
#include "AMReX_ParReduce.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/utilities/IOManager.H"
#include "amr-wind/utilities/io_utils.H"
#include "amr-wind/utilities/linear_interpolation.H"

namespace amr_wind::MetMastMarking {

namespace {} // namespace

MetMastMarking::MetMastMarking(CFDSim& sim)
    : m_sim(sim)
    , m_repo(sim.repo())
    , m_mesh(sim.mesh())
    , m_metmast_blank(sim.repo().declare_field("metmast_blank", 1, 1, 1))
{
    amrex::ParmParse pp_abl("ABL");
    pp_abl.query("metmast_location_file", m_metmast_file);
    if (!m_metmast_file.empty()) {
        std::ifstream mastfile(m_metmast_file, std::ios::in);
        if (!mastfile.good()) {
            amrex::Abort("Cannot find Met Mast profile file " + m_metmast_file);
        }
        //! x y z horizontal_radius vertical_radius
        amrex::Real value1, value2, value3, value4, value5;
        while (mastfile >> value1 >> value2 >> value3 >> value4 >> value5) {
            m_metmast_x.push_back(value1);
            m_metmast_y.push_back(value2);
            m_metmast_z.push_back(value3);
            m_metmast_horizontal_radius.push_back(value4);
            m_metmast_vertical_radius.push_back(value5);
        }
    } else {
        amrex::Abort("Cannot find Met Mast profile file " + m_metmast_file);
    }
    int num_wind_values = static_cast<int>(m_metmast_x.size());
    m_metmast_x_d.resize(num_wind_values);
    m_metmast_y_d.resize(num_wind_values);
    m_metmast_z_d.resize(num_wind_values);
    m_metmast_horizontal_radius_d.resize(num_wind_values);
    m_metmast_vertical_radius_d.resize(num_wind_values);
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_metmast_x.begin(), m_metmast_x.end(),
        m_metmast_x_d.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_metmast_y.begin(), m_metmast_y.end(),
        m_metmast_y_d.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_metmast_z.begin(), m_metmast_z.end(),
        m_metmast_z_d.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_metmast_horizontal_radius.begin(),
        m_metmast_horizontal_radius.end(),
        m_metmast_horizontal_radius_d.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_metmast_vertical_radius.begin(),
        m_metmast_vertical_radius.end(), m_metmast_vertical_radius_d.begin());
    m_sim.io_manager().register_output_var("metmast_blank");
    m_metmast_blank.setVal(0.0);
    m_metmast_blank.set_default_fillpatch_bc(m_sim.time());
}

void MetMastMarking::initialize_fields(int level, const amrex::Geometry& geom)
{

    BL_PROFILE("amr-wind::" + this->identifier() + "::initialize_fields");

    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();
    auto& blanking = m_metmast_blank(level);
    auto levelBlanking = blanking.arrays();
    const auto* metmast_x_d = m_metmast_x_d.data();
    const auto* metmast_y_d = m_metmast_y_d.data();
    const auto* metmast_z_d = m_metmast_z_d.data();
    const auto* horizontal_radius_d = m_metmast_horizontal_radius_d.data();
    const auto* vertical_radius_d = m_metmast_vertical_radius_d.data();
    for (int ii = 0; ii < static_cast<int>(m_metmast_x.size()); ii++) {
        amrex::Print() << "Met Mast " << ii << " at " << metmast_x_d[ii] << ", "
                       << metmast_y_d[ii] << ", " << metmast_z_d[ii]
                       << " with horizontal radius " << horizontal_radius_d[ii]
                       << " and vertical radius " << vertical_radius_d[ii]
                       << "\n";
        amrex::ParallelFor(
            blanking, m_metmast_blank.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
                const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
                const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
                amrex::Real ri2 = (x - metmast_x_d[ii]) * (x - metmast_x_d[ii]);
                ri2 += (y - metmast_y_d[ii]) * (y - metmast_y_d[ii]);
                if (((ri2 <=
                      (horizontal_radius_d[ii] * horizontal_radius_d[ii])) &&
                     (std::abs(z - metmast_z_d[ii]) <=
                      vertical_radius_d[ii])) ||
                    (levelBlanking[nbx](i, j, k, 0) == 1.0)) {
                    levelBlanking[nbx](i, j, k, 0) = 1.0;
                } else {
                    levelBlanking[nbx](i, j, k, 0) = 0.0;
                }
            });
    }
    amrex::Gpu::streamSynchronize();
}

void MetMastMarking::post_regrid_actions()
{
    const int nlevels = m_sim.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        initialize_fields(lev, m_sim.repo().mesh().Geom(lev));
    }
}

} // namespace amr_wind::MetMastMarking

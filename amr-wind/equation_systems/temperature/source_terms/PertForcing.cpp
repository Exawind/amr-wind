

#include "amr-wind/equation_systems/temperature/source_terms/PertForcing.H"
#include "amr-wind/utilities/IOManager.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include "AMReX_Random.H"

namespace amr_wind::pde::temperature {

PertForcing::PertForcing(const CFDSim& sim)
    : m_sim(sim), m_time(sim.time()), m_mesh(sim.mesh())
{
    amrex::ParmParse pp("PertForcing");
    pp.query("xstart", m_xstart);
    pp.query("xend", m_xend);
    pp.query("ystart", m_ystart);
    pp.query("yend", m_yend);
    pp.query("zstart", m_zstart);
    pp.query("zend", m_zend);
    pp.query("time_steps", m_time_index);
    pp.query("start_level", m_start_level);
    pp.query("pert_amplitude", m_pert_amplitude);
}

PertForcing::~PertForcing() = default;

void PertForcing::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState /*fstate */,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto& index = m_time.time_index();
    const bool has_terrain =
        this->m_sim.repo().int_field_exists("terrain_blank");
    const auto* m_terrain_blank =
        has_terrain ? &this->m_sim.repo().get_int_field("terrain_blank")
                    : nullptr;
    const auto& blank_arr = has_terrain
                                ? (*m_terrain_blank)(lev).const_array(mfi)
                                : amrex::Array4<int>();
    const auto* m_terrain_height =
        has_terrain ? &this->m_sim.repo().get_field("terrain_height") : nullptr;
    const auto& height_arr = has_terrain
                                 ? (*m_terrain_height)(lev).const_array(mfi)
                                 : amrex::Array4<double>();
    if (index % m_time_index == 0 && index > 0 && lev > m_start_lev) {
        const auto& geom = m_mesh.Geom(lev);
        const auto& dx = geom.CellSizeArray();
        const auto& prob_lo = geom.ProbLoArray();
        const amrex::Real xstart = m_xstart;
        const amrex::Real ystart = m_ystart;
        const amrex::Real zstart = m_zstart;
        const amrex::Real xend = m_xend;
        const amrex::Real yend = m_yend;
        const amrex::Real zend = m_zend;
        amrex::Real pert_amplitude = static_cast<float>(
            2.0 * m_pert_amplitude * (amrex::Random() - 0.5));
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            const amrex::Real height_arr_cell =
                (has_terrain) ? height_arr(i, j, k, 0) : 0.0;
            const amrex::Real blank_arr_cell =
                (has_terrain) ? blank_arr(i, j, k, 0) : 0.0;
            const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
            const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
            const amrex::Real z = std::max(
                prob_lo[2] + (k + 0.5) * dx[2] - height_arr_cell, 0.5 * dx[2]);
            if (x >= xstart && x <= xend && y >= ystart && y <= yend &&
                z >= zstart && z <= zend) {
                src_term(i, j, k, 0) += (1.0 - blank_arr_cell) * pert_amplitude;
            }
        });
    }
}

} // namespace amr_wind::pde::temperature

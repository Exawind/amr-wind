#include "amr-wind/equation_systems/temperature/source_terms/DragTempForcing.H"
#include "amr-wind/utilities/IOManager.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include "AMReX_Random.H"

namespace amr_wind::pde::temperature {

DragTempForcing::DragTempForcing(const CFDSim& sim)
    : m_mesh(sim.mesh())
    , m_temperature(sim.repo().get_field("temperature"))
    , m_terrainBlank(sim.repo().get_field("terrainBlank"))
{
    amrex::ParmParse pp("DragTempForcing");
    pp.query("dragCoefficient", m_drag);
}

DragTempForcing::~DragTempForcing() = default;

void DragTempForcing::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState fstate,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto temperature = m_temperature(lev).const_array(mfi);
    const auto blank = m_terrainBlank(lev).const_array(mfi);
    const auto& geom_vec = m_mesh.Geom();
    const auto& geom = geom_vec[lev];
    const auto& dx = geom.CellSize();
    const amrex::Real gpu_drag = m_drag;
    const amrex::Real gpu_TRef = 300.0;
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        amrex::Real Cd = gpu_drag / dx[0];
        src_term(i, j, k, 0) -=
            (Cd * (temperature(i, j, k, 0) - gpu_TRef) * blank(i, j, k));
    });
}

} // namespace amr_wind::pde::temperature

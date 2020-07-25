#include "amr-wind/overset/TiogaInterface.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/equation_systems/PDEBase.H"

namespace amr_wind {

namespace {

/** Convert iblanks to AMReX mask
 *
 *  \f{align}
 *  \mathrm{mask}_{i,j,k} = \begin{cases}
 *  0 & \mathrm{IBLANK}_{i, j, k} > 0 \\
 *  1 & \mathrm{IBLANK}_{i, j, k} \leq 0
 *  \end{cases}
 *  \f}
 */
void iblank_to_mask(const IntField& iblank, IntField& maskf)
{
    const auto& nlevels = iblank.repo().mesh().finestLevel() + 1;

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& ibl = iblank(lev);
        auto& mask = maskf(lev);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(ibl); mfi.isValid(); ++mfi) {
            const auto& gbx = mfi.growntilebox();
            const auto& ibarr = ibl.const_array(mfi);
            const auto& marr = mask.array(mfi);
            amrex::ParallelFor(
                gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    marr(i, j, k) = ibarr(i, j, k) > 0 ? 0 : 1;
                });
        }
    }
}
} // namespace

// clang-format off

TiogaInterface::TiogaInterface(CFDSim& sim)
    : m_iblank_cell(sim.repo().declare_int_field(
          "iblank_cell", 1, sim.pde_manager().num_ghost_state()))
    , m_iblank_node(sim.repo().declare_int_field(
          "iblank_node", 1, sim.pde_manager().num_ghost_state(), 1,
          FieldLoc::NODE))
    , m_mask_cell(sim.repo().declare_int_field(
          "mask_cell", 1, sim.pde_manager().num_ghost_state()))
    , m_mask_node(sim.repo().declare_int_field(
          "mask_node", 1, sim.pde_manager().num_ghost_state(), 1,
          FieldLoc::NODE))
{}
// clang-format on

void TiogaInterface::post_init_actions() {}

void TiogaInterface::post_regrid_actions() {}

void TiogaInterface::pre_overset_conn_work()
{
    m_iblank_cell.setVal(1);
    m_iblank_node.setVal(1);
}

void TiogaInterface::post_overset_conn_work()
{
    iblank_to_mask(m_iblank_cell, m_mask_cell);
    iblank_to_mask(m_iblank_node, m_mask_node);
}

void TiogaInterface::register_solution() {}

void TiogaInterface::update_solution() {}

} // namespace amr_wind

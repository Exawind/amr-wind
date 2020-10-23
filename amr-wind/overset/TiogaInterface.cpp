#include "amr-wind/overset/TiogaInterface.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/equation_systems/PDEBase.H"
#include "amr-wind/core/field_ops.H"
#include "amr-wind/utilities/IOManager.H"

namespace amr_wind {

namespace {

/** Convert iblanks to AMReX mask
 *
 *  \f{align}
 *  \mathrm{mask}_{i,j,k} = \begin{cases}
 *  1 & \mathrm{IBLANK}_{i, j, k} = 0 \\
 *  0 & \mathrm{IBLANK}_{i, j, k} \leq 0
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
                    marr(i, j, k) = amrex::max(ibarr(i, j, k), 0);
                });
        }
    }
}
} // namespace

// clang-format off

TiogaInterface::TiogaInterface(CFDSim& sim)
    : m_sim(sim)
    , m_iblank_cell(sim.repo().declare_int_field(
          "iblank_cell", 1, sim.pde_manager().num_ghost_state()))
    , m_iblank_node(sim.repo().declare_int_field(
          "iblank_node", 1, sim.pde_manager().num_ghost_state(), 1,
          FieldLoc::NODE))
    , m_mask_cell(sim.repo().declare_int_field(
          "mask_cell", 1, sim.pde_manager().num_ghost_state()))
    , m_mask_node(sim.repo().declare_int_field(
          "mask_node", 1, sim.pde_manager().num_ghost_state(), 1,
          FieldLoc::NODE))
{
    m_sim.io_manager().register_output_int_var(m_iblank_cell.name());
}
// clang-format on

void TiogaInterface::post_init_actions()
{
    amr_to_tioga_mesh();

    // Initialize masking so that all cells are active in solvers
    m_mask_cell.setVal(1);
    m_mask_node.setVal(1);
}

void TiogaInterface::post_regrid_actions()
{
    amr_to_tioga_mesh();

    // Initialize masking so that all cells are active in solvers
    m_mask_cell.setVal(1);
    m_mask_node.setVal(1);
}

void TiogaInterface::pre_overset_conn_work()
{
    m_iblank_cell.setVal(1);
    m_iblank_node.setVal(1);
}

void TiogaInterface::post_overset_conn_work()
{
    iblank_to_mask(m_iblank_cell, m_mask_cell);
    iblank_to_mask(m_iblank_node, m_mask_node);

    // Update equation systems after a connectivity update
    m_sim.pde_manager().icns().post_regrid_actions();
    for (auto& eqn : m_sim.pde_manager().scalar_eqns())
        eqn->post_regrid_actions();
}

void TiogaInterface::register_solution()
{
    m_qcell.reset();
    m_qnode.reset();
    auto& repo = m_sim.repo();
    auto& vel = repo.get_field("velocity");
    auto& pres = repo.get_field("p");
    m_qcell = repo.create_scratch_field(
        vel.num_comp(), vel.num_grow()[0], vel.field_location());
    m_qnode = repo.create_scratch_field(
        pres.num_comp(), pres.num_grow()[0], pres.field_location());

    vel.fillpatch(0.0);
    pres.fillpatch(0.0);

    field_ops::copy(*m_qcell, vel, 0, 0, vel.num_comp(), vel.num_grow());
    field_ops::copy(*m_qnode, pres, 0, 0, pres.num_comp(), pres.num_grow());
}

void TiogaInterface::update_solution()
{
    auto& repo = m_sim.repo();
    auto& vel = repo.get_field("velocity").state(amr_wind::FieldState::Old);
    auto& pres = repo.get_field("p");

    field_ops::copy(vel, *m_qcell, 0, 0, vel.num_comp(), vel.num_grow());
    field_ops::copy(pres, *m_qnode, 0, 0, pres.num_comp(), pres.num_grow());

#if 0
    int nlevels = repo.num_active_levels();
    for (int lev=0; lev < nlevels; ++lev) {
        auto& vfab = vel(lev);
        auto& pfab = pres(lev);
        auto& qcfab = (*m_qcell)(lev);
        auto& qnfab = (*m_qnode)(lev);
        for (amrex::MFIter mfi(vfab); mfi.isValid(); ++mfi) {
            {
                const int ncomp = vel.num_comp();
                auto bx = mfi.tilebox();
                auto varr = vfab.array(mfi);
                auto qcarr = qcfab.const_array(mfi);
                auto marr = m_mask_cell(lev).const_array(mfi);

                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        amrex::Real fac =
                            static_cast<amrex::Real>(marr(i, j, k));
                        for (int n = 0; n < ncomp; ++n) {
                            varr(i, j, k, n) = fac * varr(i, j, k, n) +
                                               (1.0 - fac) * qcarr(i, j, k, n);
                        }
                    });
            }
            {
                auto bx = mfi.nodaltilebox();
                auto parr = pfab.array(mfi);
                auto qnarr = qnfab.const_array(mfi);
                auto marr = m_mask_node(lev).const_array(mfi);

                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        amrex::Real fac =
                            static_cast<amrex::Real>(marr(i, j, k));
                        parr(i, j, k) =
                            fac * parr(i, j, k) + (1.0 - fac) * qnarr(i, j, k);
                    });
            }
        }
    }
#endif

    // fixme this is only necessary because tioga does
    // not fill in ghosts yet
    vel.fillpatch(0.0);
    pres.fillpatch(0.0);
}

void TiogaInterface::amr_to_tioga_mesh()
{
    BL_PROFILE("amr-wind::TiogaInterface::amr_to_tioga_mesh");
    auto& mesh = m_sim.mesh();
    const int nlevels = mesh.finestLevel() + 1;
    const int iproc = amrex::ParallelDescriptor::MyProc();
    const int nproc = amrex::ParallelDescriptor::NProcs();
    const auto* problo = mesh.Geom(0).ProbLo();

    int ngrids_global = 0;
    int ngrids_local = 0;
    for (int lev = 0; lev < nlevels; ++lev) {
        ngrids_global += mesh.boxArray(lev).size();

        const auto& dmap = mesh.DistributionMap(lev);
        for (long d = 0; d < dmap.size(); ++d) {
            if (dmap[d] == iproc) ++ngrids_local;
        }
    }

    // Reset array data structures
    m_info.ngrids_global = ngrids_global;
    m_info.ngrids_local = ngrids_local;
    m_info.int_data.resize(AMROversetInfo::ints_per_grid * ngrids_global);
    m_info.real_data.resize(AMROversetInfo::reals_per_grid * ngrids_global);
    m_info.gid_map.resize(nlevels);
    for (auto& vv : m_info.gid_map) vv.clear();

    // Track local grid ID
    std::vector<int> lgrid_id(nproc, 0);

    auto& idata = m_info.int_data;
    auto& rdata = m_info.real_data;

    int igp = 0; // Global index of the grid
    int iix = 0; // Index into the integer array
    int irx = 0; // Index into the real array

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& ba = mesh.boxArray(lev);
        const auto& dm = mesh.DistributionMap(lev);
        const amrex::Real* dx = mesh.Geom(lev).CellSize();

        for (long d = 0; d < dm.size(); ++d) {
            idata[iix++] = igp;             // Global index of this patch
            idata[iix++] = lev;             // Level of this patch
            idata[iix++] = dm[d];           // MPI rank of this patch
            idata[iix++] = lgrid_id[dm[d]]; // Local ID for this patch

            const auto& bx = ba[d];
            const int* lo = bx.loVect();
            const int* hi = bx.hiVect();

            for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                idata[iix + i] = lo[i];
                idata[iix + AMREX_SPACEDIM + i] = hi[i];

                rdata[irx + i] = problo[i] + lo[i] * dx[i];
                rdata[irx + AMREX_SPACEDIM + i] = dx[i];
            }
            iix += 2 * AMREX_SPACEDIM;
            irx += 2 * AMREX_SPACEDIM;

            if (iproc == dm[d]) {
                m_info.gid_map[lev].push_back(igp);
            }

            // Increment global ID counter
            ++igp;
            // Increment local index
            ++lgrid_id[dm[d]];
        }
    }
}

} // namespace amr_wind

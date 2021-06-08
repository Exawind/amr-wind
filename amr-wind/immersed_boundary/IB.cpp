#include "amr-wind/immersed_boundary/IB.H"
#include "amr-wind/immersed_boundary/IBModel.H"
#include "amr-wind/immersed_boundary/IBParser.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/FieldRepo.H"

#include <algorithm>

namespace amr_wind {
namespace ib {

IB::IB(CFDSim& sim) : m_sim(sim) {}

IB::~IB() = default;

void IB::pre_init_actions()
{
    BL_PROFILE("amr-wind::ib::IB::post_init_actions");
    amrex::ParmParse pp(identifier());

    amrex::Vector<std::string> labels;
    pp.getarr("labels", labels);

    const int n_ibs = labels.size();

    for (int i = 0; i < n_ibs; ++i) {
        const std::string& tname = labels[i];
        const std::string& prefix = identifier() + "." + tname;
        amrex::ParmParse pp1(prefix);

        std::string type;
        pp.query("type", type);
        pp1.query("type", type);
        AMREX_ALWAYS_ASSERT(!type.empty());

        auto obj = ImmersedBoundaryModel::create(type, m_sim, tname, i);

        const std::string default_prefix = identifier() + "." + type;
        utils::IBParser inp(default_prefix, prefix);

        obj->read_inputs(inp);
        m_ibs.emplace_back(std::move(obj));
    }
}

void IB::post_init_actions()
{
    BL_PROFILE("amr-wind::ib::IB::post_init_actions");
    for (auto& ib : m_ibs) ib->init_ib();
}

void IB::post_regrid_actions() {}

void IB::pre_advance_work()
{
    BL_PROFILE("amr-wind::ib::IB::pre_advance_work");
}

void IB::pre_pressure_correction_work()
{
    BL_PROFILE("amr-wind::ib::IB::pre_pressure_correction_work");
    update_velocities();
}

/** Update immersed boundary positions.
 *
 *  This method loops over all immersed boundaries and updates
 *  their position. Their new location may be prescribed or computed
 *
 *  \sa IB::update_velocities
 */
void IB::update_positions()
{
    BL_PROFILE("amr-wind::ib::IB::update_positions");
}

/** Provide updated velocities to immersed boundary instances
 *
 *  \sa IB::update_positions
 */
void IB::update_velocities()
{
    BL_PROFILE("amr-wind::ib::IB::update_velocities");
    const int nlevels = m_sim.repo().num_active_levels();
    auto& mask_cell = m_sim.repo().get_int_field("mask_cell");
    auto& velocity = m_sim.repo().get_field("velocity");

    for (int lev = 0; lev < nlevels; ++lev) {
        for (amrex::MFIter mfi(mask_cell(lev)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.growntilebox();
            auto epsilon_cell = mask_cell(lev).array(mfi);
            auto varr = velocity(lev).array(mfi);
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    if (epsilon_cell(i, j, k) == 0) {
                        varr(i, j, k, 0) = 0.;
                        varr(i, j, k, 1) = 0.;
                        varr(i, j, k, 2) = 0.;
                    }
                });
        }
    }
}

/** Helper method to compute forces on all immersed boundary components
 */
void IB::compute_forces()
{
    BL_PROFILE("amr-wind::ib::IB::compute_forces");
    for (auto& ib : m_ibs) {
        ib->compute_forces();
    }
}

void IB::prepare_outputs()
{
    const std::string out_dir_prefix = "post_processing/immersed_boundary";
    const std::string sname =
        amrex::Concatenate(out_dir_prefix, m_sim.time().time_index());
    if (!amrex::UtilCreateDirectory(sname, 0755)) {
        amrex::CreateDirectoryFailed(sname);
    }
    for (auto& ib : m_ibs) {
        ib->prepare_outputs(sname);
    }
}

void IB::post_advance_work()
{
    for (auto& ib : m_ibs) {
        ib->write_outputs();
    }
}

} // namespace ib
} // namespace amr_wind

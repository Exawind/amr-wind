#include "amr-wind/immersed_boundary/IB.H"
#include "amr-wind/immersed_boundary/IBModel.H"
#include "amr-wind/immersed_boundary/IBParser.H"
#include "amr-wind/immersed_boundary/IBContainer.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/FieldRepo.H"

#include <algorithm>

namespace amr_wind {
namespace ib {

IB::IB(CFDSim& sim)
    : m_sim(sim), m_ib_source(sim.repo().declare_field("ib_src_term", 3))
{}

IB::~IB() = default;

void IB::pre_init_actions()
{
    BL_PROFILE("amr-wind::ib::IB::post_init_actions");
    amrex::ParmParse pp(identifier());

    amrex::Vector<std::string> labels;
    pp.getarr("labels", labels);

    const int nturbines = labels.size();

    for (int i = 0; i < nturbines; ++i) {
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

    amrex::Vector<int> ib_proc_count(amrex::ParallelDescriptor::NProcs(), 0);
    for (auto& ib : m_ibs) ib->determine_root_proc(ib_proc_count);

    {
        // Sanity check that we have processed the turbines correctly
        int nib =
            std::accumulate(ib_proc_count.begin(), ib_proc_count.end(), 0);
        AMREX_ALWAYS_ASSERT(num_ibs() == nib);
    }

    for (auto& ib : m_ibs) ib->init_ib_source();

    setup_container();
    update_positions();
    update_velocities();
    compute_forces();
    compute_source_term();
    prepare_outputs();
}

void IB::post_regrid_actions()
{
    for (auto& ib : m_ibs) ib->determine_influenced_procs();

    setup_container();
}

void IB::pre_advance_work()
{
    BL_PROFILE("amr-wind::ib::IB::pre_advance_work");

    m_container->reset_container();
    update_positions();
    update_velocities();
    compute_forces();
    compute_source_term();
}

/** Set up the container for sampling velocities
 *
 *  Allocates memory and initializes the particles corresponding to immersed boundary
 *  nodes for all turbines that influence the current MPI rank. This method is
 *  invoked once during initialization and during regrid step.
 */
void IB::setup_container()
{
    const int ntotal = num_ibs();
    const int nlocal = std::count_if(
        m_ibs.begin(), m_ibs.end(),
        [](const std::unique_ptr<ImmersedBoundaryModel>& obj) {
            return obj->info().sample_vel_in_proc;
        });

    m_container.reset(new IBContainer(m_sim.mesh(), nlocal));

    auto& pinfo = m_container->m_data;
    for (int i = 0, il = 0; i < ntotal; ++i) {
        if (m_ibs[i]->info().sample_vel_in_proc) {
            pinfo.global_id[il] = i;
            pinfo.num_pts[il] = m_ibs[i]->num_velocity_points();
            ++il;
        }
    }

    m_container->initialize_container();
}

/** Update immersed boundary positions and sample velocities at new locations.
 *
 *  This method loops over all the turbines local to this MPI rank and updates
 *  the position vectors. These new locations are provided to the sampling
 *  container that samples velocities at these new locations.
 *
 *  \sa IB::update_velocities
 */
void IB::update_positions()
{
    BL_PROFILE("amr-wind::ib::IB::update_positions");
    auto& pinfo = m_container->m_data;
    for (int i = 0, ic = 0; i < pinfo.num_objects; ++i) {
        const auto ig = pinfo.global_id[i];
        auto vpos =
            ::amr_wind::utils::slice(pinfo.position, ic, pinfo.num_pts[i]);
        m_ibs[ig]->update_positions(vpos);
        ic += pinfo.num_pts[i];
    }
    m_container->update_positions();

    // Sample velocities at the new locations
    auto& vel = m_sim.repo().get_field("velocity");
    m_container->sample_velocities(vel);
}

/** Provide updated velocities from container to immersed boundary instances
 *
 *  \sa IB::update_positions
 */
void IB::update_velocities()
{
    BL_PROFILE("amr-wind::ib::IB::update_velocities");
    auto& pinfo = m_container->m_data;
    for (int i = 0, ic = 0; i < pinfo.num_objects; ++i) {
        const auto ig = pinfo.global_id[i];
        const auto vel =
            ::amr_wind::utils::slice(pinfo.velocity, ic, pinfo.num_pts[i]);
        m_ibs[ig]->update_velocities(vel);
        ic += pinfo.num_pts[i];
    }
}

/** Helper method to compute forces on all immersed boundary components
 */
void IB::compute_forces()
{
    BL_PROFILE("amr-wind::ib::IB::compute_forces");
    for (auto& ac : m_ibs) {
        if (ac->info().ib_in_proc) {
            ac->compute_forces();
        }
    }
}

void IB::compute_source_term()
{
    BL_PROFILE("amr-wind::ib::IB::compute_source_term");
    m_ib_source.setVal(0.0);
    const int nlevels = m_sim.repo().num_active_levels();

    for (int lev = 0; lev < nlevels; ++lev) {
        auto& sfab = m_ib_source(lev);
        const auto& geom = m_sim.mesh().Geom(lev);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(sfab); mfi.isValid(); ++mfi) {
            for (auto& ac : m_ibs) {
                if (ac->info().ib_in_proc) {
                    ac->compute_source_term(lev, mfi, geom);
                }
            }
        }
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
    const int iproc = amrex::ParallelDescriptor::MyProc();
    for (auto& ac : m_ibs) {
        if (ac->info().root_proc == iproc) {
            ac->prepare_outputs(sname);
        }
    }
}

void IB::post_advance_work()
{
    const int iproc = amrex::ParallelDescriptor::MyProc();
    for (auto& ac : m_ibs) {
        if (ac->info().root_proc == iproc) {
            ac->write_outputs();
        }
    }
}

} // namespace ib
} // namespace amr_wind

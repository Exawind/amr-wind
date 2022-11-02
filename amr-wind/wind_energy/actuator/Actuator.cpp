#include "amr-wind/wind_energy/actuator/Actuator.H"
#include "amr-wind/wind_energy/actuator/ActuatorModel.H"
#include "amr-wind/wind_energy/actuator/ActParser.H"
#include "amr-wind/wind_energy/actuator/ActuatorContainer.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/FieldRepo.H"

#include <algorithm>
#include <memory>

namespace amr_wind::actuator {

Actuator::Actuator(CFDSim& sim)
    : m_sim(sim), m_act_source(sim.repo().declare_field("actuator_src_term", 3))
{}

Actuator::~Actuator() = default;

void Actuator::pre_init_actions()
{
    BL_PROFILE("amr-wind::actuator::Actuator::post_init_actions");
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

        auto obj = ActuatorModel::create(type, m_sim, tname, i);

        const std::string default_prefix = identifier() + "." + type;
        utils::ActParser inp(default_prefix, prefix);

        obj->read_inputs(inp);
        m_actuators.emplace_back(std::move(obj));
    }
}

void Actuator::post_init_actions()
{
    BL_PROFILE("amr-wind::actuator::Actuator::post_init_actions");

    amrex::Vector<int> act_proc_count(amrex::ParallelDescriptor::NProcs(), 0);
    for (auto& act : m_actuators) {
        act->determine_root_proc(act_proc_count);
    }

    {
        // Sanity check that we have processed the turbines correctly
        int nact =
            std::accumulate(act_proc_count.begin(), act_proc_count.end(), 0);
        AMREX_ALWAYS_ASSERT(num_actuators() == nact);
    }

    for (auto& act : m_actuators) {
        act->init_actuator_source();
    }

    setup_container();
    update_positions();
    update_velocities();
    compute_forces();
    compute_source_term();
    prepare_outputs();
}

void Actuator::post_regrid_actions()
{
    for (auto& act : m_actuators) {
        act->determine_influenced_procs();
    }

    setup_container();
}

void Actuator::pre_advance_work()
{
    BL_PROFILE("amr-wind::actuator::Actuator::pre_advance_work");

    m_container->reset_container();
    update_positions();
    update_velocities();
    compute_forces();
    compute_source_term();
}

/** Set up the container for sampling velocities
 *
 *  Allocates memory and initializes the particles corresponding to actuator
 *  nodes for all turbines that influence the current MPI rank. This method is
 *  invoked once during initialization and during regrid step.
 */
void Actuator::setup_container()
{
    const int ntotal = num_actuators();
    const int nlocal = std::count_if(
        m_actuators.begin(), m_actuators.end(),
        [](const std::unique_ptr<ActuatorModel>& obj) {
            return obj->info().sample_vel_in_proc;
        });

    m_container = std::make_unique<ActuatorContainer>(m_sim.mesh(), nlocal);

    auto& pinfo = m_container->m_data;
    for (int i = 0, il = 0; i < ntotal; ++i) {
        if (m_actuators[i]->info().sample_vel_in_proc) {
            pinfo.global_id[il] = i;
            pinfo.num_pts[il] = m_actuators[i]->num_velocity_points();
            ++il;
        }
    }

    m_container->initialize_container();
}

/** Update actuator positions and sample velocities at new locations.
 *
 *  This method loops over all the turbines local to this MPI rank and updates
 *  the position vectors. These new locations are provided to the sampling
 *  container that samples velocities at these new locations.
 *
 *  \sa Actuator::update_velocities
 */
void Actuator::update_positions()
{
    BL_PROFILE("amr-wind::actuator::Actuator::update_positions");
    auto& pinfo = m_container->m_data;
    for (int i = 0, ic = 0; i < pinfo.num_objects; ++i) {
        const auto ig = pinfo.global_id[i];
        auto vpos =
            ::amr_wind::utils::slice(pinfo.position, ic, pinfo.num_pts[i]);
        m_actuators[ig]->update_positions(vpos);
        ic += pinfo.num_pts[i];
    }
    m_container->update_positions();

    // Sample velocities at the new locations
    const auto& vel = m_sim.repo().get_field("velocity");
    const auto& density = m_sim.repo().get_field("density");
    m_container->sample_fields(vel, density);
}

/** Provide updated velocities from container to actuator instances
 *
 *  \sa Acuator::update_positions
 */
void Actuator::update_velocities()
{
    BL_PROFILE("amr-wind::actuator::Actuator::update_velocities");
    auto& pinfo = m_container->m_data;
    for (int i = 0, ic = 0; i < pinfo.num_objects; ++i) {
        const auto ig = pinfo.global_id[i];

        const auto vel =
            ::amr_wind::utils::slice(pinfo.velocity, ic, pinfo.num_pts[i]);

        const auto density =
            ::amr_wind::utils::slice(pinfo.density, ic, pinfo.num_pts[i]);

        m_actuators[ig]->update_fields(vel, density);
        ic += pinfo.num_pts[i];
    }
}

/** Helper method to compute forces on all actuator components
 */
void Actuator::compute_forces()
{
    BL_PROFILE("amr-wind::actuator::Actuator::compute_forces");
    for (auto& ac : m_actuators) {
        if (ac->info().actuator_in_proc) {
            ac->compute_forces();
        }
    }
}

void Actuator::compute_source_term()
{
    BL_PROFILE("amr-wind::actuator::Actuator::compute_source_term");
    m_act_source.setVal(0.0);
    const int nlevels = m_sim.repo().num_active_levels();

    for (int lev = 0; lev < nlevels; ++lev) {
        auto& sfab = m_act_source(lev);
        const auto& geom = m_sim.mesh().Geom(lev);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(sfab); mfi.isValid(); ++mfi) {
            for (auto& ac : m_actuators) {
                if (ac->info().actuator_in_proc) {
                    ac->compute_source_term(lev, mfi, geom);
                }
            }
        }
    }
}

void Actuator::prepare_outputs()
{
    const std::string out_dir_prefix = "post_processing/actuator";
    const std::string sname =
        amrex::Concatenate(out_dir_prefix, m_sim.time().time_index());
    if (!amrex::UtilCreateDirectory(sname, 0755)) {
        amrex::CreateDirectoryFailed(sname);
    }
    const int iproc = amrex::ParallelDescriptor::MyProc();
    for (auto& ac : m_actuators) {
        if (ac->info().root_proc == iproc) {
            ac->prepare_outputs(sname);
        }
    }
}

void Actuator::post_advance_work()
{
    const int iproc = amrex::ParallelDescriptor::MyProc();
    for (auto& ac : m_actuators) {
        if (ac->info().root_proc == iproc) {
            ac->write_outputs();
        }
    }
}

} // namespace amr_wind::actuator

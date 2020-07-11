#ifndef CFDSIM_H
#define CFDSIM_H

#include "AMReX_AmrCore.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/equation_systems/PDEBase.H"
#include "amr-wind/core/Physics.H"

/** AMR-Wind
 *
 *  All C++ code in AMR-Wind is organized within the amr_wind namespace.
 */
namespace amr_wind {
class IOManager;
class PostProcessManager;

namespace turbulence {
class TurbulenceModel;
}

/**
 *  \defgroup core AMR-Wind core
 *  Core data structures in AMR-Wind
 *
 *  This group describes the core data structures implemented in AMR-Wind that
 *  provide a higher-level abstraction over the classes available in AMReX.
 */

/** Data structures for a CFD simulation
 *  \ingroup core
 *
 *  CFDSim represents the amr-wind execution environment and manages all the
 *  necessary components used to perform a simulation. Each simulation contains
 *  a single CFDSim instance that holds references to the mesh, FieldRepo,
 *  SimTime, PhysicsMgr, pde::PDEMgr, IOManager, and post-processing manager
 *  instances. This class is just a data holder and does not perform any
 *  computational logic itself. The computational heavy-lifting is performed by
 *  the class instances within a time integration loop.
 */
class CFDSim
{
public:
    CFDSim(amrex::AmrCore& mesh);

    ~CFDSim();

    //! Return the AMR mesh hierarchy
    amrex::AmrCore& mesh() { return m_mesh; }
    const amrex::AmrCore& mesh() const { return m_mesh; }

    //! Return simulation time control
    SimTime& time() { return m_time; }
    const SimTime& time() const { return m_time; }

    //! Return the field repository
    FieldRepo& repo() { return m_repo; }
    FieldRepo& repo() const { return m_repo; }

    pde::PDEMgr& pde_manager() { return m_pde_mgr; }
    const pde::PDEMgr& pde_manager() const { return m_pde_mgr; }

    PhysicsMgr& physics_manager() { return m_physics_mgr; }
    const PhysicsMgr& physics_manager() const { return m_physics_mgr; }

    PhysicsMgr::TypeVector& physics() { return m_physics_mgr.objects(); }
    const PhysicsMgr::TypeVector& physics() const { return m_physics_mgr.objects(); }

    turbulence::TurbulenceModel& turbulence_model() { return *m_turbulence; }
    const turbulence::TurbulenceModel& turbulence_model() const
    { return *m_turbulence; }

    IOManager& io_manager() { return *m_io_mgr; }
    const IOManager& io_manager() const { return * m_io_mgr; }

    PostProcessManager& post_manager() { return *m_post_mgr; }
    const PostProcessManager& post_manager() const  { return *m_post_mgr; }

    //! Instantiate the turbulence model based on user inputs
    void create_turbulence_model();

    //! Instantiate all physics instances based on user inputs
    void init_physics();

private:
    amrex::AmrCore& m_mesh;

    SimTime m_time;

    mutable FieldRepo m_repo;

    pde::PDEMgr m_pde_mgr;

    PhysicsMgr m_physics_mgr;

    std::unique_ptr<turbulence::TurbulenceModel> m_turbulence;

    std::unique_ptr<IOManager> m_io_mgr;

    std::unique_ptr<PostProcessManager> m_post_mgr;
};

}

#endif /* CFDSIM_H */

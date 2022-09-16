#include "amr-wind/CFDSim.H"
#include "amr-wind/turbulence/TurbulenceModel.H"
#include "amr-wind/utilities/IOManager.H"
#include "amr-wind/utilities/PostProcessing.H"
#include "amr-wind/overset/OversetManager.H"
#include "amr-wind/core/ExtSolver.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {

CFDSim::CFDSim(amrex::AmrCore& mesh)
    : m_mesh(mesh)
    , m_repo(m_mesh)
    , m_pde_mgr(*this)
    , m_io_mgr(new IOManager(*this))
    , m_post_mgr(new PostProcessManager(*this))
    , m_ext_solver_mgr(new ExtSolverMgr)
    , m_helics(new helics_storage(*this)) {}

CFDSim::~CFDSim() = default;

void CFDSim::create_turbulence_model()
{
    std::string transport_model = "ConstTransport";
    std::string turbulence_model = "Laminar";
    {
        amrex::ParmParse pp("transport");
        pp.query("model", transport_model);
    }
    {
        amrex::ParmParse pp("turbulence");
        pp.query("model", turbulence_model);
    }

    const std::string identifier = turbulence_model + "-" + transport_model;
    m_turbulence = turbulence::TurbulenceModel::create(identifier, *this);
    m_turbulence->parse_model_coeffs();
}

void CFDSim::init_physics()
{
    amrex::ParmParse pp("incflo");
    amrex::Vector<std::string> phys_names;
    pp.queryarr("physics", phys_names);

    for (auto& phy : phys_names) {
        m_physics_mgr.create(phy, *this);
    }
}

void CFDSim::activate_overset()
{
    amrex::ParmParse pp("overset");
    std::string otype = "TIOGA";

    m_overset_mgr = OversetManager::create(otype, *this);
}

bool CFDSim::has_overset() const { return (static_cast<bool>(m_overset_mgr)); }

void CFDSim::activate_mesh_map()
{
    amrex::ParmParse pp("geometry");
    std::string mesh_map_name; // default
    m_mesh_mapping = static_cast<bool>(pp.query("mesh_mapping", mesh_map_name));
    if (m_mesh_mapping) {
        m_mesh_map = MeshMap::create(mesh_map_name);
        m_mesh_map->declare_mapping_fields(*this, m_pde_mgr.num_ghost_state());
    }
}

} // namespace amr_wind

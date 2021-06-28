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
{}

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
    // declare nodal and cell-centered mesh mapping array
    auto& mesh_scale_fac_cc = m_repo.declare_cc_field(
        "mesh_scaling_factor_cc", AMREX_SPACEDIM, 0, 1);
    auto& mesh_scale_fac_nd = m_repo.declare_nd_field(
        "mesh_scaling_factor_nd", AMREX_SPACEDIM, 0, 1);

    amrex::ParmParse pp("geometry");
    if(pp.contains("mesh_mapping")) {
        std::string mesh_map_name;
        pp.query("mesh_mapping", mesh_map_name);
        m_mesh_map_mgr.create(mesh_map_name, *this);
    }
    else {
        // always create default mesh mapping
        m_mesh_map_mgr.create("ConstantScaling", *this);
    }
}

} // namespace amr_wind

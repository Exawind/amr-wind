#include "amr-wind/CFDSim.H"
#include "amr-wind/turbulence/TurbulenceModel.H"
#include "amr-wind/utilities/IOManager.H"
#include "amr-wind/utilities/PostProcessing.H"
#include "amr-wind/overset/OversetManager.H"

#ifdef AMREX_USE_DPCPP
#include "amr-wind/turbulence/turbulence_utils.H"
#endif

#include "AMReX_ParmParse.H"

namespace amr_wind {

CFDSim::CFDSim(amrex::AmrCore& mesh)
    : m_mesh(mesh)
    , m_time()
    , m_repo(m_mesh)
    , m_pde_mgr(*this)
    , m_io_mgr(new IOManager(*this))
    , m_post_mgr(new PostProcessManager(*this))
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

#ifndef AMREX_USE_DPCPP
    const std::string identifier = turbulence_model + "-" + transport_model;
    m_turbulence = turbulence::TurbulenceModel::create(identifier, *this);
#else
    m_turbulence = turbulence::create_turbulence_model(turbulence_model, *this);
#endif
}

void CFDSim::init_physics()
{
    amrex::ParmParse pp("incflo");
    amrex::Vector<std::string> phys_names;
    pp.queryarr("physics", phys_names);

    for (auto& phy: phys_names)
        m_physics_mgr.create(phy, *this);
}

void CFDSim::activate_overset()
{
    amrex::ParmParse pp("overset");
    std::string otype = "TIOGA";

    m_overset_mgr = OversetManager::create(otype, *this);
}

} // namespace amr_wind

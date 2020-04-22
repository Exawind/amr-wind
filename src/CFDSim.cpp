#include "CFDSim.H"
#include "TurbulenceModel.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {

CFDSim::CFDSim(amrex::AmrCore& mesh)
    : m_mesh(mesh), m_time(), m_repo(m_mesh), m_pde_mgr(*this)
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
}

} // namespace amr_wind

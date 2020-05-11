#include "PDEBase.H"
#include "SchemeTraits.H"
#include "CFDSim.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace pde {

PDEMgr::PDEMgr(CFDSim& sim)
    : m_sim(sim)
    , m_probtype(0)
{
    amrex::ParmParse pp("incflo");
    pp.query("probtype", m_probtype);
    pp.query("use_godunov", m_use_godunov);

    m_scheme = m_use_godunov
        ? fvm::Godunov::scheme_name()
        : fvm::MOL::scheme_name();
}

PDEBase& PDEMgr::register_icns()
{
    const std::string name = "ICNS-" + m_scheme;

    m_icns = amr_wind::pde::PDEBase::create(name, m_sim, m_probtype);
    return *m_icns;
}

PDEBase& PDEMgr::register_transport_pde(const std::string& pde_name)
{
    const std::string name = pde_name + "-" + m_scheme;

    if (contains(name)) {
        amrex::Print()
            << "WARNING: Multiple requests to register PDE: " << pde_name
            << std::endl;
        return operator()(name);
    }

    return create(name, m_sim, m_probtype);
}

bool PDEMgr::has_pde(const std::string &pde_name) const
{
    const std::string name = pde_name + "-" + m_scheme;
    return contains(name);
}

}
}

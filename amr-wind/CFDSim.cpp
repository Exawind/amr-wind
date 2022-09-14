#include "amr-wind/CFDSim.H"
#include "amr-wind/turbulence/TurbulenceModel.H"
#include "amr-wind/utilities/IOManager.H"
#include "amr-wind/utilities/PostProcessing.H"
#include "amr-wind/overset/OversetManager.H"
#include "amr-wind/core/ExtSolver.H"

#include "AMReX_ParmParse.H"


#include "helics/cpp98/CombinationFederate.hpp"
#include "helics/cpp98/helics.hpp"
#include "helics/cpp98/Federate.hpp"


#include <sstream>
#include <iostream>
#include <string>

using namespace helicscpp;

namespace amr_wind {

helics_storage::helics_storage()
	: m_fi(new helicscpp::FederateInfo("zmq"))
    , m_vfed(new helicscpp::CombinationFederate("Test receiver Federate",*m_fi))
{

    std::cout <<"PI RECEIVER: Value federate created";
    m_vfed->setProperty(HELICS_PROPERTY_TIME_DELTA, 1.0);
    //  Subscribe to PI SENDER's publication 
    // auto& in1 = vFed1->registerInput<double>("");
    std::cout <<"PI RECEIVER: Value federate created";
    // auto& subid = cFed->registerSubscription(vtarget + "/pub", "double");
    auto m_sub_2 = m_vfed->registerSubscription("control", "string");
    // m_sub = m_vfed->registerSubscription("control", "string");
    //std::cout <<"PI RECEIVER: Subscription registered";
    std::cout <<"PI RECEIVER: Value federate created";
    //  Register the publication 
    auto m_pub_2 = m_vfed->registerGlobalPublication("status", "string");


    //std::cout <<"AMRWIND RECEIVER: Publication registered\n";
    std::cout <<"PI RECEIVER: Value federate created";
    //  Enter initialization state 
    m_vfed->enterInitializingMode(); // can throw helicscpp::InvalidStateTransition exception
    std::cout <<"PI RECEIVER: Value federate created";
    m_vfed->enterExecutingMode(); 
    std::cout <<"PI RECEIVER: Value federate created";
    // currenttime = m_vfed->getCurrentTime();
    std::cout <<"Creation complete!! ";
    // std::cout <<"\n" << m_sub.getString();
    
    m_inflow_wind_speed_to_amrwind = 0.0;
    m_inflow_wind_direction_to_amrwind = 0.0;    
    m_turbine_power_to_controller.resize(m_num_turbines, 0.0);
    m_turbine_yaw_to_controller.resize(m_num_turbines, 0.0);  
    m_turbine_yaw_to_amrwind.resize(m_num_turbines, 0.0);  
        
}

helics_storage::~helics_storage() = default;

CFDSim::CFDSim(amrex::AmrCore& mesh)
    : m_mesh(mesh)
    , m_repo(m_mesh)
    , m_pde_mgr(*this)
    , m_io_mgr(new IOManager(*this))
    , m_post_mgr(new PostProcessManager(*this))
    , m_ext_solver_mgr(new ExtSolverMgr)
    , m_helics(new helics_storage) {}

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

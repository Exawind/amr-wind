#include "amr-wind/helics.H"

#include <sstream>
#include <iostream>
#include <string>

using namespace helicscpp;

namespace amr_wind {

helics_storage::helics_storage(CFDSim& sim)
	: m_fi(new helicscpp::FederateInfo("zmq"))
    , m_vfed(new helicscpp::CombinationFederate("Test receiver Federate",*m_fi))
    , m_sim(sim)
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

void helics_storage::send_messages_to_controller() 
{

    // send time, speed (1 + 1)
    // send turbine powers to controller (num_trubines)
    if (amrex::ParallelDescriptor::IOProcessor())
    {
    	// put helics send stuff here
    }
}

void helics_storage::recv_messages_from_controller() 
{

    // receive wind direction and speed from controller (1 + 1)
	// receive turbine yaw directions (num_turbines)
	if (amrex::ParallelDescriptor::IOProcessor())
    {
        std::stringstream charFromControlCenter;
        HelicsTime currenttime = m_vfed->requestNextStep();
        std::cout << "\n error at 94";
        int subCount = m_vfed->getInputCount();
        std::cout << "\n error at 102   "<<subCount;
        helicscpp::Input sub;

        for(int i = 0; i < subCount; i++) {
            int yy = i;
            sub = m_vfed->getSubscription(yy);
            std::cout << "\n advancing time "<< sub.getString().c_str();
        }

        std::stringstream ssToControlCenter; 

        helicscpp::Publication pub; 

        int pubCount = m_vfed->getPublicationCount();

        for(int i = 0; i < pubCount; i++) {
            int yy = i;
            pub = m_vfed->getPublication(yy);
            ssToControlCenter << "[all random values form amrwind!! "<<currenttime<<"]";
            std::string strToControlCenter; 
            strToControlCenter = ssToControlCenter.str();
            pub.publish(strToControlCenter.c_str());
        }
    }
    
    
	// broadcast wind turbine yaw angles to all procs
	// FIXME: some day only need to send/recv to specific turbines
    amrex::ParallelDescriptor::Bcast(
            m_turbine_yaw_to_amrwind.data(), m_turbine_yaw_to_amrwind.size(),
            amrex::ParallelDescriptor::IOProcessorNumber(),
            amrex::ParallelDescriptor::Communicator());
        
    // broadcast wind speed and direction to all procs
    amrex::ParallelDescriptor::Bcast(
            &m_inflow_wind_speed_to_amrwind, 1,
            amrex::ParallelDescriptor::IOProcessorNumber(),
            amrex::ParallelDescriptor::Communicator());
            
    amrex::ParallelDescriptor::Bcast(
            &m_inflow_wind_direction_to_amrwind, 1,
            amrex::ParallelDescriptor::IOProcessorNumber(),
            amrex::ParallelDescriptor::Communicator());

}

helics_storage::~helics_storage() = default;

}

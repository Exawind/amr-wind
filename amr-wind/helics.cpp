#include "amr-wind/helics.H"

#include "AMReX_ParmParse.H"

#include <sstream>
#include <iostream>
#include <string>
#include <stdlib.h>

#include <list>

using namespace helicscpp;

void tokenize(std::string s, std::string del, std::list<double> &return_list)
{

  int start = 0;
  int end = s.find(del);

  return_list.push_front(atof(s.substr(start + 1, end - start).c_str()));
  while (end > 0)
  {

    start = end + del.size();
    end = s.find(del, start);

    if (end == -1)
      end = -2;
    return_list.push_front(atof(s.substr(start, end - start).c_str()));
  }
}



namespace amr_wind {

helics_storage::helics_storage(CFDSim& sim)
    : m_sim(sim)
{

    if (amrex::ParallelDescriptor::IOProcessor()) {
	
		m_fi = std::make_unique<helicscpp::FederateInfo>("zmq");
		m_vfed = std::make_unique<helicscpp::CombinationFederate>("Test receiver Federate",*m_fi);    
	
        m_vfed->setProperty(HELICS_PROPERTY_TIME_DELTA, 1.0);
        //  Subscribe to PI SENDER's publication
        m_vfed->registerSubscription("control", "string");
        //  Register the publication
        m_vfed->registerGlobalPublication("status", "string");
        //  Enter initialization state
        m_vfed->enterInitializingMode(); // can throw
                                         // helicscpp::InvalidStateTransition
                                         // exception
        m_vfed->enterExecutingMode();
        std::cout << "PI RECEIVER: Value federate created";
        std::cout << "Creation complete!! ";
    }
    
    // parse input file and count how many actuators exist
    amrex::Vector<std::string> actuators;
    amrex::ParmParse pp("Actuator");
    pp.queryarr("labels", actuators);

    if (actuators.size() > 0) {
        m_num_turbines = actuators.size();
    }

    m_turbine_power_to_controller.resize(m_num_turbines, 0.0);
    m_turbine_yaw_to_controller.resize(m_num_turbines, 0.0);
    m_turbine_yaw_to_amrwind.resize(m_num_turbines, 0.0);
}

void helics_storage::send_messages_to_controller()
{

    // send time, speed (1 + 1)
    // send turbine powers to controller (num_trubines)
    if (amrex::ParallelDescriptor::IOProcessor()) {
        // put helics send stuff here
    }
}

void helics_storage::recv_messages_from_controller()
{

    // receive wind direction and speed from controller (1 + 1)
    // receive turbine yaw directions (num_turbines)
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::stringstream charFromControlCenter;
        currenttime = m_vfed->requestNextStep();
        int subCount = m_vfed->getInputCount();
        helicscpp::Input sub;

        for (int i = 0; i < subCount; i++) {
            int yy = i;
            sub = m_vfed->getSubscription(yy);
            charFromControlCenter << sub.getString().c_str();
            std::cout << "\n advancing time " << charFromControlCenter.str();
        }


        if (currenttime > 1)
        {
            // Igonre timestep 0 since message pipe has junk. 
            // // Unpack the values from the control center using a string stream

			std::list<double> return_list;
			tokenize(charFromControlCenter.str(), ",", return_list);

			m_inflow_wind_direction_to_amrwind = return_list.front();
			return_list.pop_front();
			m_inflow_wind_speed_to_amrwind = return_list.front();
			return_list.pop_front();
			float time= return_list.front();
			return_list.pop_front();

			std::cout << "\n speed: "<<m_inflow_wind_speed_to_amrwind <<"  direction: "<<m_inflow_wind_direction_to_amrwind<<" Time "<<time << std::endl;

        }

        std::stringstream ssToControlCenter;

        helicscpp::Publication pub;

        int pubCount = m_vfed->getPublicationCount();

        for (int i = 0; i < pubCount; i++) {
            pub = m_vfed->getPublication(i);

            ssToControlCenter << "[" ;

            for (int yy=0;yy<m_turbine_power_to_controller.size(); yy++)
            {    
            	ssToControlCenter <<m_turbine_power_to_controller[yy];
            	ssToControlCenter << ",";
            }
            
            ssToControlCenter << m_sim.time().current_time() << "]";

            std::string strToControlCenter;
            strToControlCenter = ssToControlCenter.str();
            pub.publish(strToControlCenter.c_str());
            
            std::cout<<"\n should send m_turbine_power_to_controller "<<strToControlCenter.c_str()<<std::endl;


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


} // namespace amr_wind

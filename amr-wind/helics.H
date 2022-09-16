#ifndef HELICS_H
#define HELICS_H

#include "AMReX_AmrCore.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/CFDSim.H"

#include "helics/cpp98/CombinationFederate.hpp"
#include "helics/cpp98/helics.hpp"
#include "helics/cpp98/Federate.hpp"


namespace amr_wind {

class helics_storage
{
public: 
	
	explicit helics_storage(CFDSim& sim);
	~helics_storage();
	
public:

	void send_messages_to_controller();
	void recv_messages_from_controller();
	
    HelicsTime currenttime = 0.0;
    std::string fedinitstring = "--federates=1";
    double deltat = 0.01; // what is this for? should we use m_time.deltaT() instead?

    std::unique_ptr<helicscpp::FederateInfo> m_fi;
    std::unique_ptr<helicscpp::CombinationFederate> m_vfed;
    
    const int m_num_turbines = 100;
    
    double m_inflow_wind_speed_to_amrwind;
    double m_inflow_wind_direction_to_amrwind;
    
    std::vector<double> m_turbine_power_to_controller;
    std::vector<double> m_turbine_yaw_to_controller;    
    std::vector<double> m_turbine_yaw_to_amrwind;    
    
private:

    CFDSim& m_sim;
    
};


} // namespace amr_wind

#endif /* HELICS_H */
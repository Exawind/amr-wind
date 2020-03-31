#include "abl_test_utils.H"

namespace amr_wind_tests {
namespace utils {

void populate_abl_params()
{
    amrex::ParmParse pp("abl");

    // Initial conditions (Temperature)
    amrex::Vector<amrex::Real> theights{{0.0, 650.0, 750.0, 1000.0}};
    amrex::Vector<amrex::Real> tvalues{{300.0, 300.0, 308.0, 308.75}};
    pp.add("ntemperature", theights.size());
    pp.addarr("temperature_heights", theights);
    pp.addarr("temperature_values", tvalues);
    pp.add("perturb_ref_height", 50.0);

    // Boussinesq Buoyancy
    pp.add("reference_temperature", 300.0);

    // ABL Forcing
    pp.add("abl_forcing_height", 90.0);

    // Coriolis term
    pp.add("latitude", 45.0);

    pp.add("kappa", 0.41);
    pp.add("surface_roughness_z0", 0.1);

    // Needed for initial conditions
    {
        amrex::ParmParse pp("incflo");
        pp.add("ro_0", 1.0);   // Density
        pp.add("ic_u", 20.0);  // Ux
        pp.add("ic_v", 10.0);   // Uy
        pp.add("ic_w", 0.0);   // Uw
        amrex::Vector<amrex::Real> grav{{0.0,0.0,-9.81}};
        pp.addarr("gravity",grav);
        
    }

    // Adjust computational domain to be more like ABL mesh in the z direction
    {
        amrex::ParmParse pp("amr");
        amrex::Vector<int> ncell{{8, 8, 64}};
        pp.addarr("n_cell", ncell);
    }
    {
        amrex::ParmParse pp("geometry");
        amrex::Vector<amrex::Real> probhi {{120.0, 120.0, 1000.0}};

        pp.addarr("prob_hi", probhi);
    }
}

}
}

#include <incflo.H> 

void incflo::ReadRheologyParameters()
{
     amrex::ParmParse pp("incflo");

     fluid_model = "newtonian";
     pp.query("fluid_model", fluid_model);
     if(fluid_model == "newtonian")
     {
         amrex::Print() << "Newtonian fluid with"
                        << " mu = " << mu << std::endl;
     }
     else if(fluid_model == "powerlaw")
     {
         pp.query("n", n_0);
         AMREX_ALWAYS_ASSERT(n_0 > 0.0);
         AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_0 != 1.0,
                 "No point in using power-law rheology with n = 1");

         amrex::Print() << "Power-law fluid with"
                        << " mu = " << mu
                        << ", n = " << n_0 <<  std::endl;
     }
     else if(fluid_model == "bingham")
     {
         pp.query("tau_0", tau_0);
         AMREX_ALWAYS_ASSERT_WITH_MESSAGE(tau_0 > 0.0,
                 "No point in using Bingham rheology with tau_0 = 0");

         pp.query("papa_reg", papa_reg);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(papa_reg > 0.0,
                    "Papanastasiou regularisation parameter must be positive");

         amrex::Print() << "Bingham fluid with"
                        << " mu = " << mu
                        << ", tau_0 = " << tau_0
                        << ", papa_reg = " << papa_reg << std::endl;
     }
     else if(fluid_model == "hb")
     {
         pp.query("n", n_0);
         AMREX_ALWAYS_ASSERT(n_0 > 0.0);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_0 != 1.0,
                 "No point in using Herschel-Bulkley rheology with n = 1");

         pp.query("tau_0", tau_0);
         AMREX_ALWAYS_ASSERT_WITH_MESSAGE(tau_0 > 0.0,
                 "No point in using Herschel-Bulkley rheology with tau_0 = 0");

         pp.query("papa_reg", papa_reg);
         AMREX_ALWAYS_ASSERT_WITH_MESSAGE(papa_reg > 0.0,
                 "Papanastasiou regularisation parameter must be positive");

         amrex::Print() << "Herschel-Bulkley fluid with"
                        << " mu = " << mu
                        << ", n = " << n_0
                        << ", tau_0 = " << tau_0
                        << ", papa_reg = " << papa_reg << std::endl;
     }
     else if(fluid_model == "smd")
     {
         pp.query("n", n_0);
         AMREX_ALWAYS_ASSERT(n_0 > 0.0);

         pp.query("tau_0", tau_0);
         AMREX_ALWAYS_ASSERT_WITH_MESSAGE(tau_0 > 0.0,
                 "No point in using de Souza Mendes-Dutra rheology with tau_0 = 0");

         pp.query("eta_0", eta_0);
         AMREX_ALWAYS_ASSERT(eta_0 > 0.0);

         amrex::Print() << "de Souza Mendes-Dutra fluid with"
                        << " mu = " << mu
                        << ", n = " << n_0
                        << ", tau_0 = " << tau_0
                        << ", eta_0 = " << eta_0 << std::endl;
     }
     else
     {
         amrex::Abort("Unknown fluid_model! Choose either newtonian, powerlaw, bingham, hb, smd");
     }
}

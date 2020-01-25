#include <incflo.H> 

using namespace amrex;

void incflo::ReadRheologyParameters()
{
     amrex::ParmParse pp("incflo");

     std::string fluid_model_s = "newtonian";
     pp.query("fluid_model", fluid_model_s);
     if(fluid_model_s == "newtonian")
     {
         m_fluid_model = FluidModel::Newtonian;
         amrex::Print() << "Newtonian fluid with"
                        << " mu = " << m_mu << std::endl;
     }
     else if(fluid_model_s == "powerlaw")
     {
         m_fluid_model = FluidModel::powerlaw;
         pp.query("n", m_n_0);
         AMREX_ALWAYS_ASSERT(m_n_0 > 0.0);
         AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_n_0 != 1.0,
                 "No point in using power-law rheology with n = 1");

         amrex::Print() << "Power-law fluid with"
                        << " mu = " << m_mu
                        << ", n = " << m_n_0 <<  std::endl;
     }
     else if(fluid_model_s == "bingham")
     {
         m_fluid_model = FluidModel::Bingham;
         pp.query("tau_0", m_tau_0);
         AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_tau_0 > 0.0,
                 "No point in using Bingham rheology with tau_0 = 0");

         pp.query("papa_reg", m_papa_reg);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_papa_reg > 0.0,
                    "Papanastasiou regularisation parameter must be positive");

         amrex::Print() << "Bingham fluid with"
                        << " mu = " << m_mu
                        << ", tau_0 = " << m_tau_0
                        << ", papa_reg = " << m_papa_reg << std::endl;
     }
     else if(fluid_model_s == "hb")
     {
         m_fluid_model = FluidModel::HerschelBulkley;
         pp.query("n", m_n_0);
         AMREX_ALWAYS_ASSERT(m_n_0 > 0.0);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_n_0 != 1.0,
                 "No point in using Herschel-Bulkley rheology with n = 1");

         pp.query("tau_0", m_tau_0);
         AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_tau_0 > 0.0,
                 "No point in using Herschel-Bulkley rheology with tau_0 = 0");

         pp.query("papa_reg", m_papa_reg);
         AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_papa_reg > 0.0,
                 "Papanastasiou regularisation parameter must be positive");

         amrex::Print() << "Herschel-Bulkley fluid with"
                        << " mu = " << m_mu
                        << ", n = " << m_n_0
                        << ", tau_0 = " << m_tau_0
                        << ", papa_reg = " << m_papa_reg << std::endl;
     }
     else if(fluid_model_s == "smd")
     {
         m_fluid_model = FluidModel::deSouzaMendesDutra;
         pp.query("n", m_n_0);
         AMREX_ALWAYS_ASSERT(m_n_0 > 0.0);

         pp.query("tau_0", m_tau_0);
         AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_tau_0 > 0.0,
                 "No point in using de Souza Mendes-Dutra rheology with tau_0 = 0");

         pp.query("eta_0", m_eta_0);
         AMREX_ALWAYS_ASSERT(m_eta_0 > 0.0);

         amrex::Print() << "de Souza Mendes-Dutra fluid with"
                        << " mu = " << m_mu
                        << ", n = " << m_n_0
                        << ", tau_0 = " << m_tau_0
                        << ", eta_0 = " << m_eta_0 << std::endl;
     }
     else
     {
         amrex::Abort("Unknown fluid_model! Choose either newtonian, powerlaw, bingham, hb, smd");
     }
}

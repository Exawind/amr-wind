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
     else if(fluid_model_s == "SmagorinskyLillySGS")
     {
        m_fluid_model = FluidModel::SmagorinskyLillySGS;
        pp.query("SmagorinskyLillyConstant", m_Smagorinsky_Lilly_SGS_constant);
        AMREX_ALWAYS_ASSERT(m_Smagorinsky_Lilly_SGS_constant > 0.0);

        amrex::Print() << "Smagorinsky Lilly SGS model"
                       << " mu = " << m_mu 
                       << " SmagorinskyLillyConstant = " << m_Smagorinsky_Lilly_SGS_constant << std::endl;
     }
     else
     {
         amrex::Abort("Unknown fluid_model! Choose either newtonian, powerlaw, bingham, hb, smd");
     }
}

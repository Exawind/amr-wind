#include "ABL.H"
#include "ABLFieldInit.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {

ABL::ABL()
{
    amrex::ParmParse pp("abl");

    pp.query("use_boussinesq", m_has_boussinesq);
    pp.query("coriolis_effect", m_has_coriolis);
    pp.query("abl_forcing", m_has_driving_dpdx);

    // Instantiate the ABL field initializer
    m_field_init = std::make_unique<ABLFieldInit>();
}

}

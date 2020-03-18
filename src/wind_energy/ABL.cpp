#include "ABL.H"
#include "ABLFieldInit.H"

#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"

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

void ABL::initialize_fields(
    const amrex::Geometry& geom, amrex::MultiFab& density,
    amrex::MultiFab& velocity, amrex::MultiFab& /* pressure */,
    amrex::MultiFab& scalars) const
{
    for (amrex::MFIter mfi(density); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();

        (*m_field_init)(
            vbx, geom, velocity.array(mfi), density.array(mfi),
            scalars.array(mfi));
    }
}

}

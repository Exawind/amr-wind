#ifndef ABLFILLCODEDINLET_H
#define ABLFILLCODEDINLET_H

#include "amr-wind/core/FieldFillPatchOps.H"
#include "amr-wind/core/FieldBCOps.H"
#include "amr-wind/wind_energy/ABLCodedInlet.H"

namespace amr_wind {

/** Operations to perform inflow of data from boundary planes
 *  \ingroup we_abl
 *
 *  This class contains the operations necessary to interface with
 *  ABLBoundaryPlane class and populate data at the boundaries
 *
 *  \sa ABLBoundaryPlane
 */
class ABLFillCodedInlet : public FieldFillPatchOps<FieldBCDirichlet>
{
public:
    ABLFillCodedInlet(
        Field& field,
        const amrex::AmrCore& mesh,
        const SimTime& time,
        const ABLCodedInlet& abl_coded_inlet);

    ~ABLFillCodedInlet() override;

    void fillpatch(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost,
        const FieldState fstate = FieldState::New) override;

    void fillpatch_sibling_fields(
        int lev,
        amrex::Real time,
        amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& mfabs,
        amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& ffabs,
        amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& cfabs,
        const amrex::IntVect& nghost,
        const amrex::Vector<amrex::BCRec>& bcrec,
        const FieldState fstate = FieldState::New,
        const FieldInterpolator itype =
            FieldInterpolator::FaceDivFree) override;

    //! Implementation that handles filling patches from a coarse to fine level
    void fillpatch_from_coarse(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost,
        const FieldState fstate = FieldState::New) override;

    //! Implementation that handles filling physical boundary conditions
    void fillphysbc(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost,
        const FieldState fstate = FieldState::New) override;

protected:
    const ABLCodedInlet& m_abl_coded_inlet;
};

} // namespace amr_wind

#endif /* ABLFILLCODEDINLET_H */

#ifndef ABLFILLINFLOW_H
#define ABLFILLINFLOW_H

#include "amr-wind/core/FieldFillPatchOps.H"
#include "amr-wind/core/FieldBCOps.H"
#include "amr-wind/wind_energy/ABLBoundaryPlane.H"

namespace amr_wind {

/** Operations to perform inflow of data from boundary planes
 *  \ingroup we_abl
 *
 *  This class contains the operations necessary to interface with
 *  ABLBoundaryPlane class and populate data at the boundaries
 *
 *  \sa ABLBoundaryPlane
 */
class ABLFillInflow : public FieldFillPatchOps<FieldBCNoOp>
{
public:
    ABLFillInflow(
        Field& field,
        const amrex::AmrCore& mesh,
        const SimTime& time,
        const ABLBoundaryPlane& abl_bc);

    virtual ~ABLFillInflow();

    virtual void fillpatch(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost) override;

    //! Implementation that handles filling patches from a coarse to fine level
    virtual void fillpatch_from_coarse(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost) override;

    //! Implementation that handles filling physical boundary conditions
    virtual void fillphysbc(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost) override;

protected:
    const ABLBoundaryPlane& m_bndry_plane;
};

} // namespace amr_wind

#endif /* ABLFILLINFLOW_H */

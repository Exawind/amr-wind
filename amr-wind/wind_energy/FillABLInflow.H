#ifndef FILLABLINFLOW_H
#define FILLABLINFLOW_H

#include "amr-wind/core/FieldFillPatchOps.H"
#include "amr-wind/core/FieldBCOps.H"

namespace amr_wind {

class FillABLInflow : public FieldFillPatchOps<FieldBCNoOp>
{
public:
    FillABLInflow(
        Field& field, const amrex::AmrCore& mesh, const SimTime& time
        //, const ABLBoundaryPlane
    );

    virtual ~FillABLInflow();

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
    // const ABLBoundaryPlane m_abl_bc;
};

} // namespace amr_wind

#endif /* FILLABLINFLOW_H */

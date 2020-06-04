#include "amr-wind/wind_energy/FillABLInflow.H"

namespace amr_wind {

FillABLInflow::FillABLInflow(
    Field& field,
    const amrex::AmrCore& mesh,
    const SimTime& time //, const ABLBoundaryPlane
    )
    : FieldFillPatchOps<FieldBCNoOp>(
          field, mesh, time, 0, FieldInterpolator::CellConsLinear)
{}

FillABLInflow::~FillABLInflow() = default;

void FillABLInflow::fillpatch(
    int lev,
    amrex::Real time,
    amrex::MultiFab& mfab,
    const amrex::IntVect& nghost)
{
    FieldFillPatchOps<FieldBCNoOp>::fillpatch(lev, time, mfab, nghost);

    // apply ABL inflow
}

void FillABLInflow::fillpatch_from_coarse(
    int lev,
    amrex::Real time,
    amrex::MultiFab& mfab,
    const amrex::IntVect& nghost)
{
    FieldFillPatchOps<FieldBCNoOp>::fillpatch_from_coarse(
        lev, time, mfab, nghost);

    // apply ABL inflow
}

void FillABLInflow::fillphysbc(
    int lev,
    amrex::Real time,
    amrex::MultiFab& mfab,
    const amrex::IntVect& nghost)
{
    FieldFillPatchOps<FieldBCNoOp>::fillphysbc(lev, time, mfab, nghost);

    // apply ABL inflow
}

} // namespace amr_wind

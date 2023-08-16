#include "amr-wind/wind_energy/ABLFillInflow.H"

namespace amr_wind {

ABLFillInflow::ABLFillInflow(
    Field& field,
    const amrex::AmrCore& mesh,
    const SimTime& time,
    const ABLBoundaryPlane& bndry_plane)
    : FieldFillPatchOps<FieldBCDirichlet>(
          field, mesh, time, FieldInterpolator::CellConsLinear)
    , m_bndry_plane(bndry_plane)
{}

ABLFillInflow::~ABLFillInflow() = default;

void ABLFillInflow::fillpatch(
    int lev,
    amrex::Real time,
    amrex::MultiFab& mfab,
    const amrex::IntVect& nghost,
    const FieldState fstate)
{
    FieldFillPatchOps<FieldBCDirichlet>::fillpatch(
        lev, time, mfab, nghost, fstate);

    m_bndry_plane.populate_data(lev, time, m_field, mfab);
}

void ABLFillInflow::fillpatch_from_coarse(
    int lev,
    amrex::Real time,
    amrex::MultiFab& mfab,
    const amrex::IntVect& nghost,
    const FieldState fstate)
{
    FieldFillPatchOps<FieldBCDirichlet>::fillpatch_from_coarse(
        lev, time, mfab, nghost, fstate);

    m_bndry_plane.populate_data(lev, time, m_field, mfab);
}

void ABLFillInflow::fillphysbc(
    int lev,
    amrex::Real time,
    amrex::MultiFab& mfab,
    const amrex::IntVect& nghost,
    const FieldState fstate)
{
    FieldFillPatchOps<FieldBCDirichlet>::fillphysbc(
        lev, time, mfab, nghost, fstate);

    m_bndry_plane.populate_data(lev, time, m_field, mfab);
}

void ABLFillInflow::fillpatch_sibling_fields(
    int lev,
    amrex::Real time,
    amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& mfabs,
    amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& ffabs,
    amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& cfabs,
    const amrex::IntVect& nghost,
    const amrex::Vector<amrex::BCRec>& bcrec,
    const FieldState fstate,
    const FieldInterpolator itype)
{
    // For an ABL fill, we just foextrap the mac velocities
    amrex::Vector<amrex::BCRec> lbcrec(m_field.num_comp());
    const auto& ibctype = m_field.bc_type();
    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        const auto side = ori.faceDir();
        const auto bct = ibctype[ori];
        const int dir = ori.coordDir();
        for (int i = 0; i < m_field.num_comp(); ++i) {
            if (bct == BC::mass_inflow) {
                if (side == amrex::Orientation::low) {
                    lbcrec[i].setLo(dir, amrex::BCType::foextrap);
                } else {
                    lbcrec[i].setHi(dir, amrex::BCType::foextrap);
                }
            } else {
                if (side == amrex::Orientation::low) {
                    lbcrec[i].setLo(dir, bcrec[i].lo(dir));
                } else {
                    lbcrec[i].setHi(dir, bcrec[i].hi(dir));
                }
            }
        }
    }

    FieldFillPatchOps<FieldBCDirichlet>::fillpatch_sibling_fields(
        lev, time, mfabs, ffabs, cfabs, nghost, lbcrec, fstate, itype);

    for (int i = 0; i < static_cast<int>(mfabs.size()); i++) {
        m_bndry_plane.populate_data(lev, time, m_field, *mfabs[i], 0, i);
    }
}

} // namespace amr_wind

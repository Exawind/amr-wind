#include "amr-wind/wind_energy/OceanWavesFillInflow.H"

namespace amr_wind {

OceanWavesFillInflow::OceanWavesFillInflow(
    Field& field,
    const amrex::AmrCore& mesh,
    const SimTime& time)
    : FieldFillPatchOps<FieldBCDirichlet>(
          field, mesh, time, FieldInterpolator::CellConsLinear)
{}

OceanWavesFillInflow::~OceanWavesFillInflow() = default;

void OceanWavesFillInflow::fillpatch(
    int lev,
    amrex::Real time,
    amrex::MultiFab& mfab,
    const amrex::IntVect& nghost,
    const FieldState fstate)
{
    FieldFillPatchOps<FieldBCDirichlet>::fillpatch(
        lev, time, mfab, nghost, fstate);

    // ** //
}

void OceanWavesFillInflow::fillpatch_from_coarse(
    int lev,
    amrex::Real time,
    amrex::MultiFab& mfab,
    const amrex::IntVect& nghost,
    const FieldState fstate)
{
    FieldFillPatchOps<FieldBCDirichlet>::fillpatch_from_coarse(
        lev, time, mfab, nghost, fstate);

    // ** //
}

void OceanWavesFillInflow::fillphysbc(
    int lev,
    amrex::Real time,
    amrex::MultiFab& mfab,
    const amrex::IntVect& nghost,
    const FieldState fstate)
{
    FieldFillPatchOps<FieldBCDirichlet>::fillphysbc(
        lev, time, mfab, nghost, fstate);

    // ** //
}

void OceanWavesFillInflow::fillpatch_sibling_fields(
    int lev,
    amrex::Real time,
    amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& mfabs,
    amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& ffabs,
    amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& cfabs,
    const amrex::IntVect& nghost,
    const amrex::Vector<amrex::BCRec>& bcrec,
    const amrex::Vector<amrex::BCRec>& /* unused */,
    const FieldState fstate)
{
    // Avoid trying to read after planes have already been populated
    const bool plane_data_unchanged = m_bndry_plane.is_data_newer_than(time);
    // For an ABL fill, we just foextrap the mac velocities
    amrex::Vector<amrex::BCRec> fp_bcrec(m_field.num_comp());
    amrex::Vector<amrex::BCRec> ph_bcrec(m_field.num_comp());
    const auto& ibctype = m_field.bc_type();
    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        const auto side = ori.faceDir();
        const auto bct = ibctype[ori];
        const int dir = ori.coordDir();
        for (int i = 0; i < m_field.num_comp(); ++i) {
            if ((bct == BC::mass_inflow) || (bct == BC::mass_inflow_outflow)) {
                if (side == amrex::Orientation::low) {
                    ph_bcrec[i].setLo(dir, amrex::BCType::foextrap);
                    fp_bcrec[i].setLo(
                        dir, plane_data_unchanged ? amrex::BCType::ext_dir
                                                  : amrex::BCType::foextrap);
                } else {
                    ph_bcrec[i].setHi(dir, amrex::BCType::foextrap);
                    fp_bcrec[i].setHi(
                        dir, plane_data_unchanged ? amrex::BCType::ext_dir
                                                  : amrex::BCType::foextrap);
                }
            } else {
                if (side == amrex::Orientation::low) {
                    ph_bcrec[i].setLo(dir, bcrec[i].lo(dir));
                    fp_bcrec[i].setLo(dir, bcrec[i].lo(dir));
                } else {
                    ph_bcrec[i].setHi(dir, bcrec[i].hi(dir));
                    fp_bcrec[i].setHi(dir, bcrec[i].hi(dir));
                }
            }
        }
    }

    FieldFillPatchOps<FieldBCDirichlet>::fillpatch_sibling_fields(
        lev, time, mfabs, ffabs, cfabs, nghost, fp_bcrec, ph_bcrec, fstate);

    // fill target
}

} // namespace amr_wind

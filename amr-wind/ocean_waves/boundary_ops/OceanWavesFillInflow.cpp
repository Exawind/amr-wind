#include "amr-wind/ocean_waves/boundary_ops/OceanWavesFillInflow.H"

namespace amr_wind {

OceanWavesFillInflow::OceanWavesFillInflow(
    Field& field,
    const amrex::AmrCore& mesh,
    const SimTime& time,
    const OceanWavesBoundary& ow_bndry)
    : FieldFillPatchOps<FieldBCDirichlet>(
          field, mesh, time, FieldInterpolator::CellConsLinear)
    , m_ow_bndry(ow_bndry)
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

    if (m_field.base_name() == "velocity") {
        m_ow_bndry.set_velocity(lev, time, m_field, mfab);
    } else if (m_field.base_name() == "vof") {
        m_ow_bndry.set_vof(lev, time, m_field, mfab);
    } else if (m_field.base_name() == "density") {
        m_ow_bndry.set_density(lev, time, m_field, mfab);
    }
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

    if (m_field.base_name() == "velocity") {
        m_ow_bndry.set_velocity(lev, time, m_field, mfab);
    } else if (m_field.base_name() == "vof") {
        m_ow_bndry.set_vof(lev, time, m_field, mfab);
    } else if (m_field.base_name() == "density") {
        m_ow_bndry.set_density(lev, time, m_field, mfab);
    }
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

    if (m_field.base_name() == "velocity") {
        m_ow_bndry.set_velocity(lev, time, m_field, mfab);
    } else if (m_field.base_name() == "vof") {
        m_ow_bndry.set_vof(lev, time, m_field, mfab);
    } else if (m_field.base_name() == "density") {
        m_ow_bndry.set_density(lev, time, m_field, mfab);
    }
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
    // First foextrap the MAC velocities
    amrex::Vector<amrex::BCRec> lbcrec(m_field.num_comp());
    const auto& ibctype = m_field.bc_type();
    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        const auto side = ori.faceDir();
        const auto bct = ibctype[ori];
        const int dir = ori.coordDir();
        for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
            auto ori = oit();
            const auto side = ori.faceDir();
            const auto bct = ibctype[ori];
            const int dir = ori.coordDir();
            for (int i = 0; i < m_field.num_comp(); ++i) {
                if ((bct == BC::mass_inflow) ||
                    (bct == BC::mass_inflow_outflow)) {
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
    }

    FieldFillPatchOps<FieldBCDirichlet>::fillpatch_sibling_fields(
            lev, time, mfabs, ffabs, cfabs, nghost, lbcrec, lbcrec, fstate);

    for (int i = 0; i < static_cast<int>(mfabs.size()); i++) {
            m_ow_bndry.set_velocity(lev, time, m_field, *mfabs[i], 0, i);
        }
}

} // namespace amr_wind

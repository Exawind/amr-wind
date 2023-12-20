#include "amr-wind/boundary_conditions/MassOutflowBC.H"

namespace amr_wind {
namespace {
AMREX_FORCE_INLINE
amrex::Box lower_boundary_faces(const amrex::Box& b, int dir)
{
    amrex::IntVect lo(b.smallEnd());
    amrex::IntVect hi(b.bigEnd());
    int sm = lo[dir];
    lo.setVal(dir, sm - 1);
    hi.setVal(dir, sm - 1);
    amrex::IndexType bxtype(b.ixType());
    bxtype.set(dir);
    return {lo, hi, bxtype};
}
} // namespace

MassOutflowBC::MassOutflowBC(Field& field, amrex::Orientation ori)
    : m_field(field), m_ori(ori)
{}

} // namespace amr_wind
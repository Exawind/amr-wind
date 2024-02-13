
#include "amr-wind/boundary_conditions/MassInflowOutflowBC.H"

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

MassInflowOutflowBC::MassInflowOutflowBC(Field& field, amrex::Orientation ori)
    : m_field(field), m_ori(ori)
{}

void MassInflowOutflowBC::operator()(Field& /*field*/, const FieldState /*rho_state*/)
{}

} // namespace amr_wind
#include "amr-wind/physics/udfs/LinearProfile.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/core/FieldRepo.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace udf {

template <int DIM>
LinearProfile<DIM>::LinearProfile(const Field& fld) : m_op()
{
    AMREX_ALWAYS_ASSERT(DIM == fld.num_comp());
    amrex::ParmParse pp("LinearProfile." + fld.name());

    pp.query("direction", m_op.idir);
    const auto geom = fld.repo().mesh().Geom(0);
    m_op.zmin = geom.ProbLo(m_op.idir);
    m_op.zmax = geom.ProbHi(m_op.idir);
    pp.query("start", m_op.zmin);
    pp.query("stop", m_op.zmax);
    amrex::Vector<amrex::Real> start_val, end_val;
    pp.getarr("start_val", start_val);
    pp.getarr("stop_val", end_val);

    AMREX_ASSERT(start_val.size() == DIM);
    AMREX_ASSERT(end_val.size() == DIM);

    for (int i = 0; i < DIM; ++i) {
        m_op.fld_min[i] = start_val[i];
        m_op.fld_max[i] = end_val[i];
    }
}

template struct LinearProfile<1>;
template struct LinearProfile<AMREX_SPACEDIM>;

} // namespace udf
} // namespace amr_wind

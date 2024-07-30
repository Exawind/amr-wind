#include "amr-wind/physics/udfs/TwoLayer.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/core/vs/vector.H"
#include "amr-wind/equation_systems/icns/icns.H"

#include "AMReX_ParmParse.H"

namespace amr_wind::udf {

TwoLayer::TwoLayer(const Field& fld)
{
    {
       const int ncomp = fld.num_comp();
       amrex::ParmParse pp("TwoLayer");
       //pp.query("pvel", m_op.pvel);
       //pp.query("mvel", m_op.mvel);
       amrex::Vector<amrex::Real> pvel(0.0, ncomp);
       amrex::Vector<amrex::Real> mvel(0.0, ncomp);
       pp.getarr("pvel", pvel);
       pp.getarr("mvel", mvel);
       AMREX_ALWAYS_ASSERT(pvel.size() == ncomp);
       AMREX_ALWAYS_ASSERT(mvel.size() == ncomp);
       for (int i = 0; i < ncomp; ++i) {
            m_op.pvel[i] = pvel[i];
            m_op.mvel[i] = mvel[i];
       }
    }
}

} // namespace amr_wind::udf

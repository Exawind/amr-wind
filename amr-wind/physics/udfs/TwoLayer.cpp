#include "amr-wind/physics/udfs/TwoLayer.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/core/vs/vector.H"
#include "amr-wind/equation_systems/icns/icns.H"

#include "AMReX_ParmParse.H"

namespace amr_wind::udf {

TwoLayer::TwoLayer(const Field& fld)
{
    // This is a where the user can set some user defined variables
    // This capability can be activated with the following in the input file:
    // xlo.type = "mass_inflow"
    // xlo.velocity.inflow_type = CustomVelocity
    // CustomVelocity.foo = 1.0

    // clang-format off
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
    // clang-format on
    //amrex::Abort(
    //    "Please define the body of this function and the corresponding struct "
    //    "in the header file before using it. Then remove this message");
}

} // namespace amr_wind::udf

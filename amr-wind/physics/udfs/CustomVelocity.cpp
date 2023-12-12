#include "amr-wind/physics/udfs/CustomVelocity.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/core/vs/vector.H"
#include "amr-wind/equation_systems/icns/icns.H"

#include "AMReX_ParmParse.H"

namespace amr_wind::udf {

CustomVelocity::CustomVelocity(const Field& /*fld*/)
{
    // This is a where the user can set some user defined variables
    // This capability can be activated with the following in the input file:
    // xlo.type = "mass_inflow"
    // xlo.velocity.inflow_type = CustomVelocity
    // CustomVelocity.foo = 1.0

    // clang-format off
    // TODO: Do ParmParse here if needed
    //const int ncomp = fld.num_comp();
    //{
    //    amrex::ParmParse pp("incflo");
    //    amrex::Vector<amrex::Real> vel(0.0, ncomp);
    //    pp.getarr("velocity", vel);
    //    AMREX_ALWAYS_ASSERT(vel.size() == ncomp);
    //    for (int i = 0; i < ncomp; ++i) {
    //        m_op.vel_ref[i] = vel[i];
    //    }
    //}
    //{
    //    amrex::ParmParse pp("CustomVelocity");
    //    pp.query("foo", m_op.foo);
    //}
    // clang-format on
    amrex::Abort(
        "Please define the body of this function and the corresponding struct "
        "in the header file before using it. Then remove this message");
}

} // namespace amr_wind::udf
